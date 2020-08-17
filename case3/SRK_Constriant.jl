using JuMP
using MathOptInterface
const R = 8.314
const m_factor = [0.48508, 1.55171, -0.15613]
const ac_factor = 0.42747
const b_factor = 0.08664
const m = m_factor[1] .+ m_factor[2] .* ω .+ m_factor[3] .* ω.^2
const bi  = b_factor * R .* Tc ./ Pc
const ac = ac_factor * R^2 .* Tc.^2 ./ Pc




@variable(model, AL)
@variable(model, BL)
@variable(model, ZL)

@variable(model, ar[1:nc])
@variable(model, aij[i in 1:nc, j in i:nc])

@variable(model, amiL[1:nc])
@variable(model, amL)
@variable(model, bmL)
@variable(model, ϕL[i in 1:nc])


@NLconstraint(model, def_ar[i in 1:nc], ar[i] == sqrt(ac[i]) * (1 + m[i] * (1 - Temp / Tc[i])))
@NLconstraint(model, def_aij[i in 1:nc, j in i:nc], aij[i, j] == (1 - k[i,j]) * ar[i] * ar[j])
@NLconstraint(model, def_amL, amL == sum(x[i]^2 * aij[i,i] for i in 1:nc) + 2 * sum(x[i]^2 * aij[i,j] for i in 1:nc, j in i + 1:nc))
@NLconstraint(model, def_amiL[i in 1:nc], amiL[i] / x[i] == sum(aij[j, i] for j in 1:i - 1) + sum(aij[i, j] for j in i:nc))
@NLconstraint(model, def_bmL, bmL == sum(x[i] * bi[i] for i in 1:nc))

@NLconstraint(model, def_AL, AL == amL * P / R^2 / Temp^2)
@NLconstraint(model, def_BL, BL == bmL * P / R^2 / Temp^2)




@NLconstraint(model, def_phiL[i in 1:nc], ϕL[i] == (ZL - 1) * bi[i] / bmL - log(ZL - BL) - AL / BL * (2 * amiL[i] / amL  - bi[i] / bmL) * log(1 + BL / ZL))

@NLconstraint(model, def_ZL1, ZL^3 - (1 - BL) * ZL^2 + (AL - 3BL^2 - 2BL) * ZL - (AL * BL - BL^2 - BL^3) == 0)
@NLconstraint(model, def_ZL2, 3 * ZL^2 - 2 * (1 - BL) * ZL + (AL - 3BL^2 - 2BL) >= 0)
@NLconstraint(model, def_ZL3, 6 * ZL - 2 * (1 - BL)  <= 0)

@variable(model, AV)
@variable(model, BV)
@variable(model, ZV)

@variable(model, amiV[1:nc])
@variable(model, amV)
@variable(model, bmV)
@variable(model, ϕV[i in 1:nc])



@NLconstraint(model, def_AV, AV == amV * P / R^2 / Temp^2)
@NLconstraint(model, def_BV, BV == bmV * P / R^2 / Temp^2)

# @NLconstraint(model, def_arV[i in 1:nc], arV[i] == sqrt(ac[i]) * (1 + m[i] * (1 - Temp / Tc[i])))
# @NLconstraint(model, def_aij[i in 1:nc, j in i:nc], aij[i, j] == (1 - k[i,j]) * ar[i] * ar[j])
@NLconstraint(model, def_amV, amV == sum(y[i]^2 * aij[i,i] for i in 1:nc) + 2 * sum(y[i]^2 * aij[i,j] for i in 1:nc, j in i + 1:nc))
@NLconstraint(model, def_amiV[i in 1:nc], amiV[i] / y[i] == sum(aij[j, i] for j in 1:i - 1) + sum(aij[i, j] for j in i:nc))
@NLconstraint(model, def_bmV, bmV == sum(y[i] * bi[i] for i in 1:nc))

@NLconstraint(model, def_phiV[i in 1:nc], ϕV[i] == (ZV - 1) * bi[i] / bmV - log(ZV - BV) - AV / BV * (2 * amiV[i] / amV  - bi[i] / bmV) * log(1 + BV / ZV))

@NLconstraint(model, def_ZV1, ZV^3 - (1 - BV) * ZV^2 + (AV - 3BV^2 - 2BV) * ZV - (AV * BV - BV^2 - BV^3) == 0)
@NLconstraint(model, def_ZV2, 3 * ZV^2 - 2 * (1 - BV) * ZV + (AV - 3BV^2 - 2BV) >= 0)
@NLconstraint(model, def_ZV3, 6 * ZV - 2 * (1 - BV)  >= 0)

@variable(model, K[i in 1:nc])
@NLconstraint(model ,def_k[i in 1:nc], K[i] == ϕL[i] / ϕV[i])


