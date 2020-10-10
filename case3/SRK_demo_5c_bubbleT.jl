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
Vap_data = Matrix(CSV.read(joinpath(fdir, "VapPre.csv"), DataFrame)[:,2:7])
# streams = CSV.read(joinpath(fdir, "streams.csv"), DataFrame)

k_data = Matrix(df_k[2:end])       # 二元交互参数系数
Pc_data = prop.Pc .+ 101325       # 临界压力
Tc_data = prop.Tc       # 临界温度
# 偏心因子
ω_data = prop.omega
# z0_data = streams.z0
const NC_SET = collect([2:6]...)
const NC = length(NC_SET)
const VAP_PRE = Vap_data[NC_SET,:]
const K = k_data[NC_SET,NC_SET]       # 二元交互参数系数
const PRES_CRITICAL = Pc_data[NC_SET]       # 临界压力
const TEMP_CRITICAL = Tc_data[NC_SET]       # 临界温度
# 偏心因子
const ω = ω_data[NC_SET]
const PRES = 420000
const M_COFF = 10
const ρ = 500

# const PROBLEM_TYPE = "PT_FLASH"
const PROBLEM_TYPE = "BUBBLET"
const LN_EXPRESSION_A = false
const LN_EXPRESSION_PHI = false
const LC_EXPRESSION_FLAG = false
const COMPLEMENTS_FLAG = true
const USER_FUNCTION_FLAG = true
const K_MEAN_FLAG = false

const ϵx = log(10e-16)
# L0 = 0.937

z0 = [0.2, 0.2, 0.2, 0.2, 0.2]
x0 = z0
y0 = z0
F0 = 10
L0 = F0 * 0.937
V0 = F0 - L0

y0 = (F0 .* z0 - L0 .* x0) ./ V0
model = Model(optimizer_with_attributes(Ipopt.Optimizer))
include("SRK_Constriant.jl")
include("initialize.jl")

initialize(PROBLEM_TYPE;z0=z0,F0=F0,pres0=420000,temp0=56,x0=x0,y0=y0,L0=L0,V0=V0)
# initlize_bubleT(x0, PRES)
# ϕL0 =  (ZL0 - 1) .* bi ./ bmL0 .- log(ZL0 - BL0) .- AL0 ./ BL0 .* (2 .* amiL0 ./ amL0  .- bi ./ bmL0) .* log(1 + BL0 / ZL0)
# model
# @NLconstraint(model, def_y_con[i in 1:NC], y[i] == keq[i] * x[i])

# @NLobjective(model, Min, sum((y[i]  - keq[i] * x[i] )^2 for i in 1:NC))
# initialize(problem_type;z0,F0,x0,pres0,temp0)
# @NLobjective(model, Min, 100 * (sum(y[i] for i in 1:NC) - 1)^2)
# @NLconstraint(model, sumy, sum(y[i] for i in 1:NC) - 1 == 0)
# @NLobjective(model, Min, Temp)

# optimize!(model)
# println(value(Temp))
# println(sum(value.(y)))
optimize!(model)
# println.("flow L =", value(L), "  flow V=", value(V))
# # value()
# println.("sL=", value(sL), "    sV=", value(sV), "    β=", value(β))
# println(value.(ŷ + lnϕV - x̂ - lnϕL))
# value.([L * x[i] + V * y[i] - TOTAL_FLOW * z_input[NC_SET][i] for i in 1:NC])
@show SRK_EOS_Z(value.([ZL, AL, BL])...)
@show SRK_EOS_Z(value.([ZV, AV, BV])...)