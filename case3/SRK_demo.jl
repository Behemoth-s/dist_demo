using JuMP
using Ipopt
using CSV
using MathOptInterface
const MOI = MathOptInterface
fdir = @__DIR__
prop = CSV.read(joinpath(fdir, "prop.csv"))
df_k = CSV.read(joinpath(fdir, "kij.csv"))
Vp = Matrix(CSV.read(joinpath(fdir, "VapPre.csv"))[2:end])
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
const k = Matrix(df_k[2:end])
const Pc = prop.MW
const Tc = prop.Tc
const ω = prop.omega
const nc = length(Pc)

model = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1200,  "tol" => 1e-6, "print_level" => 4))
const P = 420000
using NLsolve
include("def.jl")
pres = P
temp0 = initlize_temp(pres, Vp)
K0 = initlize_K(pres, temp0, Pc, Tc, ω)

@variable(model, Temp)
@NLparameter(model, x[i in 1:nc] == x0[i])
@variable(model, y[1:nc])
include("SRK_Constriant.jl")
@NLexpression(model, def_y[i in 1:nc], y[i] - K[i] * x[i])
@NLobjective(model, Min, sum(def_y[i]^2 for i in 1:nc))

set_start_value(Temp, temp0)
set_start_value.(K, K0)
for i in 1:nc
    set_start_value(ar[i], sqrt(ac[i]) * (1 + m[i] * (1 - start_value(Temp) / Tc[i])))
end

for i in 1:nc, j in i:nc
    set_start_value(aij[i, j], (1 - k[i,j]) * start_value(ar[i]) * start_value(ar[j]))
end
let s = 0
    for i in 1:4, j in 1:4
        s += i * j
    end
end
# amL0 = sum(value(x[i])^2 * start_value(aij[i,i]) for i in 1:nc) + 2 * sum(value(x[i])^2 * start_value(aij[i, j]) for i in 1:nc, j in (i + 1):nc)

# bmL0 = sum(value(x[i]) * bi[i] for i in 1:nc)
fix(bmL, bmL0)
# @NLconstraint(model, def_aij[i in 1:nc, j in i:nc], aij[i, j] == (1 - k[i,j]) * ar[i] * ar[j])
# @NLconstraint(model, def_amL, amL == sum(x[i]^2 * aij[i,i] for i in 1:nc) + 2 * sum(x[i]^2 * aij[i,j] for i in 1:nc, j in i + 1:nc))
# @NLconstraint(model, def_amiL[i in 1:nc], amiL[i] / x[i] == sum(aij[j, i] for j in 1:i - 1) + sum(aij[i, j] for j in i:nc))
# @NLconstraint(model, def_bmL, )


# @NLconstraint(model, def_AL, AL == amL * P / R^2 / Temp^2)
# @NLconstraint(model, def_BL, BL == bmL * P / R^2 / Temp^2)




model