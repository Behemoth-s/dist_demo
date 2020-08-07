using Ipopt
using JuMP
using MathOptInterface
const MOI = MathOptInterface
using DataFrames
using CSV
nh = 22
nf = 13
nc = 5
lk = 1
hk = 5
strip_set = [2:nf - 1..., nf + 1:nh - 1...]
F0 = 1
xf0 = [0.05, 0.15, 0.25, 0.2, 0.35]
P0 = ones(nh) * 725
props = CSV.read("properties.csv")
x_ss_D = [0.11, 0.33, 0.56, 0, 0]
x_ss_W = [0, 0, 0, 0.36, 0.64]
x_ss = x_ss_D' .+ (x_ss_W - x_ss_D)' .* Vector(0:nh - 1) ./ (nh - 1)
T_ss_D = 52.37
T_ss_W = 106.11
T_ss = T_ss_D .+ (T_ss_W - T_ss_D) .* Vector(0:nh - 1) ./ (nh - 1)
Tf = 20
HF0 = sum(xf0[ j] * (props.alpha1[j] + Tf * (props.alpha2[j] + Tf * props.alpha3[j])) for j in nc)
# 创建模型
m = Model(optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 1200,  "tol" => 1e-6, "print_level" => 4))

@NLparameter(m, F == F0)
@NLparameter(m, HF == HF0)
@NLparameter(m, xf[j in 1:nc] == xf0[j])
@NLparameter(m, P[i in 1:nh] == P0[i])

@NLparameter(m, Aant[j in 1:nc] == props.Aant[j])
@NLparameter(m, Bant[j in 1:nc] == props.Bant[j])
@NLparameter(m, Cant[j in 1:nc] == props.Cant[j])

@NLparameter(m, α1[j in 1:nc] == props.alpha1[j])
@NLparameter(m, α2[j in 1:nc] == props.alpha2[j])
@NLparameter(m, α3[j in 1:nc] == props.alpha3[j])

@NLparameter(m, β1[j in 1:nc] == props.beta1[j])
@NLparameter(m, β2[j in 1:nc] == props.beta2[j])
@NLparameter(m, β3[j in 1:nc] == props.beta3[j])

@variable(m, L[1:nh] >= 0)
@variable(m, V[2:nh] >= 0)

@NLparameter(m, D == 0.42)
@NLparameter(m, RR == 1.6)

@variable(m, Qc >= 0)
@variable(m, Qr >= 0)

@variable(m, 0 <= x[1:nh, 1:nc] <= 1)
@variable(m, 0 <= y[1:nh, 1:nc] <= 1)
@variable(m, T[1:nh] >= 0)

@variable(m, HL[1:nh])
@variable(m, HV[2:nh])

# mass balance
@NLconstraint(m, mass_eq_D, L[1] + D == V[2])
@NLconstraint(m, mass_eq_strip[i in strip_set], L[i] + V[i] == L[i - 1] + V[i + 1])
@NLconstraint(m, mass_eq_W, L[nh] + V[nh] == L[nh - 1])
@NLconstraint(m, mass_eq_feed, L[nf] + V[nf] == L[nf - 1] + V[nf + 1] + F)

# energy balance
@NLconstraint(m, energy_eq_D, (D + L[1]) * HL[1] + Qc == V[2] * HV[2])
@NLconstraint(m, energy_eq_strip[i in strip_set], L[i] * HL[i] + V[i] * HV[i] == L[i - 1] * HL[i - 1] + V[i + 1] * HV[i + 1])
@NLconstraint(m, energy_eq_W, L[nh] * HL[nh] + V[nh] * HV[nh] == L[nh - 1] * HL[nh - 1] + Qr)
@NLconstraint(m, energy_eq_feed, L[nf] * HL[nf] + V[nf] * HV[nf] == L[nf - 1] * HL[nf - 1] + V[nf + 1] * HV[nf + 1] + F * HF)

# component balance
@NLconstraint(m, cmp_eq_D[j in 1:nc], x[1, j] == y[2, j])
@NLconstraint(m, cmp_eq_strip[i in strip_set, j in 1:nc], L[i] * x[i, j] + V[i] * y[i ,j] == L[i - 1] * x[i - 1, j] + V[i + 1] * y[i + 1, j]);
@NLconstraint(m, cmp_eq_W[j in 1:nc], L[nh] * x[nh, j] + V[nh] * y[nh, j] == L[nh - 1] * x[nh - 1, j])
@NLconstraint(m, cmp_eq_feed[j in 1:nc], L[nf] * x[nf, j] + V[nf] * y[nf, j] == L[nf - 1] * x[nf - 1, j] + V[nf + 1] * y[nf + 1,j] + F * xf[j]);

# phase equilibrium
@NLconstraint(m, phase_eq[i in 1:nh, j in 1:nc], y[i, j] == exp(Aant[j] + Bant[j] / (Cant[j] + T[i])) * x[i, j] / P[i]);

# normalization equation
@NLexpression(m, nm2[i in 1:nh], (sum(exp(Aant[j] + Bant[j] / (Cant[j] + T[i])) * x[i, j] for j in 1:nc) - P[i]) / P[i])

# HL
@NLconstraint(m, defHL[i in 1:nh], HL[i] == sum(x[i, j] * (α1[j] + T[i] * (α2[j] + T[i] * α3[j])) for j in nc))

# HV
@NLconstraint(m, defHV[i in 2:nh], HV[i] == sum(y[i, j] * (β1[j] + T[i] * (β2[j] + T[i] * β3[j])) for j in nc))

set_start_value.(x, x_ss)
set_start_value.(T, T_ss)

for i in 2:nh, j in 1:nc
    set_start_value(y[i,j], x_ss[i, j])
end

set_start_value.(L, 4)
set_start_value.(V, 4)

for i in 1:nh
    set_start_value(HL[i], sum(x_ss[i, j] * (props.alpha1[j] + T_ss[i] * (props.alpha2[j] + T_ss[i] * props.alpha3[j])) for j in nc))
end
for i in 2:nh
    set_start_value(HV[i], sum(x_ss[i, j] * (props.beta1[j] + T_ss[i] * (props.beta2[j] + T_ss[i] * props.beta3[j])) for j in nc))
end
@NLconstraint(m, Qclimit, Qc >= 1000)
@NLconstraint(m, defnm2[i in 1:nh], nm2[i] == 0)
@NLobjective(m, Min, (L[1] - RR * D)^2)
# @NLobjective(m, Min, Qc)
# @NLconstraint(m, rr, L[1] == RR * D)
# @NLobjective(m, Min, sum(nm2[i]^2 for i in 1:nh))
optimize!(m)
print(DataFrame(hcat(JuMP.value.(x), JuMP.value.(y), JuMP.value.(T), JuMP.value.(L)), vcat(Symbol.("x_", props[1]), Symbol.("y_", props[1]), :T, :L)))
