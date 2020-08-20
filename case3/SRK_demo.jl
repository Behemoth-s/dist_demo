using JuMP
using Ipopt
using CSV
using MathOptInterface
const MOI = MathOptInterface
using DataFrames
using NLsolve

fdir = @__DIR__

prop = CSV.read(joinpath(fdir, "prop.csv"), DataFrame)
df_k = CSV.read(joinpath(fdir, "kij.csv"), DataFrame)
const VAP_PRE = Matrix(CSV.read(joinpath(fdir, "VapPre.csv"), DataFrame)[:,2:7])
const K = Matrix(df_k[2:end])       # 二元交互参数系数
const PRES_CRITICAL = prop.Pc .+ 101325       # 临界压力
const TEMP_CRITICAL = prop.Tc       # 临界温度
# 偏心因子
const ω = prop.omega
const NC = length(PRES_CRITICAL)

const PRES = 420000

x0 = [0
0
0
4.78999999089900e-08
0.000362835999310612
0.0400588149238882
0.00975550898146453
0.00369478399297991
0.0393929919251533
0.0748969178576959
0.0460482729125083
0.00527052898998599
0.0451469959142207
0.00881364298325408
0.0372805279291670
0.0371641039293882
0.0476018539095565
0.0236012669551576
0.103495175803359
0.0102286969805655
0.0257619219510523
0.0263961829498473
0.0662105898741999
0.0333088089367133
0.0792607548494046
0.0304053059422299
0.0173125549671061
0.0114087789783233
0.00131040799751022
0.0992257598114711
0.0219731949582509
0.00331882199369424
0.0256447999512749
2.88999999450900e-07
0.0256488599512672
0
0]

model = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 6000,  "tol" => 1e-6, "print_level" => 4))
include("SRK_Constriant.jl")
include("initialize.jl")
initlize_bubleT(x0, PRES)
# ϕL0 =  (ZL0 - 1) .* bi ./ bmL0 .- log(ZL0 - BL0) .- AL0 ./ BL0 .* (2 .* amiL0 ./ amL0  .- bi ./ bmL0) .* log(1 + BL0 / ZL0)
# model
# @NLconstraint(model, def_y_con[i in 1:NC], y[i] == keq[i] * x[i])

# @NLobjective(model, Min, sum((y[i]  - keq[i] * x[i] )^2 for i in 1:NC))
@NLobjective(model, Min, 100 * (sum(y[i] for i in 1:NC) - 1)^2)
# @NLconstraint(model, sumy, sum(y[i] for i in 1:NC) - 1 == 0)
# @NLobjective(model, Min, Temp)

optimize!(model)
println(value(Temp))
println(sum(value.(y)))