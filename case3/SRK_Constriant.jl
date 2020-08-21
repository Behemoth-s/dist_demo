using JuMP
using MathOptInterface
using Complementarity
const R_VAL = 8.314e3
M_FACTOR = [0.480 1.574 0.176] # [0.48508, 1.55171, -0.15613]
const AC_FACTOR = 0.42748
const B_FACTOR = 0.08664
const MI_VAL = M_FACTOR[1] .+ M_FACTOR[2] .* ω .- M_FACTOR[3] .* ω.^2
const BI_VAL  = B_FACTOR * R_VAL .* TEMP_CRITICAL ./ PRES_CRITICAL
const AC_VAL = AC_FACTOR * R_VAL^2 .* TEMP_CRITICAL.^2 ./ PRES_CRITICAL
const u = 1
const v = 0
const I0 = sqrt(u^2 - 4v)

@variable(model, Temp)
# @variable(model, Pres)
@variable(model, x[1:NC] >= 0)
@variable(model, y[1:NC] >= 0)

@variable(model, IL1)
# @variable(model, IL2)
@variable(model, IL3)

# @variable(model, IL4)
@variable(model, IV1)

# @variable(model, IL2)
@variable(model, IV3)


@variable(model, AL)
@variable(model, BL)
@variable(model, ZL)

@variable(model, ar[1:NC])
@variable(model, aij[i in 1:NC, j in i:NC])

@variable(model, amL)
@variable(model, δL[1:NC])

@variable(model, bmL)
@variable(model, lnϕL[1:NC])

@variable(model, AV)
@variable(model, BV)
@variable(model, ZV)

@variable(model, amV)
@variable(model, δV[1:NC])

@variable(model, bmV)
@variable(model, lnϕV[1:NC])

@variable(model,sL >= 0)
@variable(model,sV >= 0)
@variable(model, β)

@NLconstraint(model, def_ar[i in 1:NC], ar[i] == sqrt(AC_VAL[i]) * (1 + MI_VAL[i] * (1 - sqrt(Temp / TEMP_CRITICAL[i]))))
@constraint(model, def_aij[i in 1:NC, j in i:NC], aij[i, j] == (1 - K[i,j]) * ar[i] * ar[j])

@NLconstraint(model, def_amL, amL == sum(x[i]^2 * aij[i,i] for i in 1:NC) + 2 * sum(x[i] * x[j] * aij[i,j] for i in 1:NC, j in i + 1:NC))
@NLconstraint(model, def_δL[i in 1:NC], δL[i] == 2 * ar[i] / amL * sum(x[j] * ar[j] * (1 - K[i,j]) for j in 1:NC))
# @NLconstraint(model, def_atmpL[i in 1:NC], atmpL[i] == sum((1 - K[i, j]) * (x[j] * ar[j]) for j in 1:NC))
# @NLconstraint(model, def_amiL[i in 1:nc], amiL[i] / x[i] == sum(aij[j, i] for j in 1:i - 1) + sum(aij[i, j] for j in i:nc))
@constraint(model, def_bmL, bmL == sum(x[i] * BI_VAL[i] for i in 1:NC))

@NLconstraint(model, def_AL, AL == amL * PRES / R_VAL^2 / Temp^2)
@NLconstraint(model, def_BL, BL == bmL * PRES / R_VAL / Temp)

# @NLconstraint(model, def_phiL[i in 1:NC], log(ϕL[i]) == (ZL - 1) * BI_VAL[i] / bmL - log(ZL - BL) - AL / BL * (2 * ar[i] * atmpL[i] / amL  - BI_VAL[i] / bmL) * log(1 + BL / ZL))

# @NLconstraint(model, def_phiL[i in 1:nc], ϕL[i] == (ZL - 1) * bi[i] / bmL - log(ZL - BL) - AL / BL * (2 * amiL[i] / amL  - bi[i] / bmL) * log(1 + BL / ZL))



# @NLconstraint(model, def_ZL3, 6 * ZL - 2   <= 0)







# @NLconstraint(model, def_arV[i in 1:nc], arV[i] == sqrt(ac[i]) * (1 + m[i] * (1 - Temp / Tc[i])))
# @NLconstraint(model, def_aij[i in 1:nc, j in i:nc], aij[i, j] == (1 - k[i,j]) * ar[i] * ar[j])
@NLconstraint(model, def_amV, amV == sum(y[i]^2 * aij[i,i] for i in 1:NC) + 2 * sum(y[i]^2 * aij[i,j] for i in 1:NC, j in i + 1:NC))
@NLconstraint(model, def_δV[i in 1:NC], δV[i] == 2 * ar[i] / amV * sum(y[j] * ar[j] * (1 - K[i,j]) for j in 1:NC))

# @NLconstraint(model, def_atmpV[i in 1:NC], atmpV[i] == sum((1 - K[i, j]) * (y[j] * ar[j]) for j in 1:NC))

# @NLconstraint(model, def_amiV[i in 1:nc], amiV[i] / y[i] == sum(aij[j, i] for j in 1:i - 1) + sum(aij[i, j] for j in i:nc))
@constraint(model, def_bmV, bmV == sum(y[i] * BI_VAL[i] for i in 1:NC))

@NLconstraint(model, def_AV, AV == amV * PRES / R_VAL^2 / Temp^2)
@NLconstraint(model, def_BV, BV == bmV * PRES / R_VAL^2 / Temp^2)

# @NLconstraint(model, def_phiV[i in 1:nc], ϕV[i] == (ZV - 1) * bi[i] / bmV - log(ZV - BV) - AV / BV * (2 * amiV[i] / amV  - bi[i] / bmV) * log(1 + BV / ZV))
# @NLconstraint(model, def_phiV[i in 1:NC], ϕV[i] == (ZV - 1) * BI_VAL[i] / bmV - log(ZV - BV) - AV / BV * (2 * ar[i] / amV * atmpV[i] - BI_VAL[i] / bmV) * log(1 + BL / ZL))
# @NLconstraint(model, def_phiV[i in 1:NC], log(ϕV[i]) == (ZV - 1) * BI_VAL[i] / bmV - log(ZV - BV) - AV / BV * (2 * ar[i] * atmpV[i] / amV  - BI_VAL[i] / bmV) * log(1 + BL / ZL))

# 简化lnϕ
@NLconstraint(model,def_IL1, exp(IL1) == ZL - BL)
@NLconstraint(model, def_IL3, exp(IL3) == ( ZL + BL ) / (ZL))
@NLconstraint(model,def_IV1, exp(IV1) == ZV - BV)
@NLconstraint(model, def_IV3, exp(IV3) == ( ZV + BV ) / (ZV))
@NLconstraint(model, def_phiL[i in 1:NC], lnϕL[i] == (ZL - 1) * BI_VAL[i] / bmL - IL1 - AL / BL * (δL[i]  - BI_VAL[i] / bmL) * IL3)
@NLconstraint(model, def_phiV[i in 1:NC], lnϕV[i] == (ZV - 1) * BI_VAL[i] / bmV - IV1 - AV / BV * (δV[i]  - BI_VAL[i] / bmV) * IV3)

# 确定三次方程根
@NLconstraint(model, def_ZL1, ZL^3 - ZL^2 + (AL - BL^2 - BL) * ZL - AL * BL == 0)
@NLconstraint(model, def_ZL2, 3 * ZL^2 - 2 *  ZL + (AL - BL^2 - BL) >= 0)
@NLconstraint(model, def_ZV1, ZV^3 - ZV^2 + (AV - BV^2 - BV) * ZV - AV * BV == 0)
@NLconstraint(model, def_ZV2, 3 * ZV^2 - 2 * ZV + (AV - BV^2 - BV) >= 0)
@NLconstraint(model, def_ZV3, 6 * ZV - 2   >= -M_COFF * sV)
@NLconstraint(model, def_ZL3, 6 * ZL - 2   <= M_COFF * sL)

# @NLconstraint(model, def_ZV3, 6 * ZV - 2   >= 0)
# @NLconstraint(model, def_sumy, sum(y[i] for i in 1:NC) == 1)





# @NLconstraint(model, def_y[i in 1:NC], y[i] == β * ϕL[i] / ϕV[i] * x[i])
@variable(model, L >= 0)
@variable(model, V >= 0)
@variable(model, x̂[1:NC])
@variable(model, ŷ[1:NC])



@NLconstraint(model, def_lnx[i in 1:NC], x[i] == exp(x̂[i]))
@NLconstraint(model, def_lny[i in 1:NC], y[i] == exp(ŷ[i]))
@constraint(model, sumy, sum(y[i] - x[i] for i in 1:NC) == 0 )

@constraint(model, lower_limit_β, β - 1  >= -sL)
@constraint(model, upper_limit_β, β - 1  <= sV)
# @constraint(model, lower_limit_β, β - 1 == sV - sL)
@NLconstraint(model, phase_eq[i in 1:NC], ŷ[i] + lnϕV[i] == log(β) + lnϕL[i] + x̂[i])

@variable(model, lc[i in 1:NC] >= 0)
@variable(model, vc[i in 1:NC] >= 0)
@constraint(model, def_lc[i in 1:NC], lc[i] == L * x[i])
@constraint(model, def_vc[i in 1:NC], vc[i] == V * y[i])

@NLexpression(model, ζ, sum(lc[i] * (lnϕL[i] + x̂[i]) + vc[i] * (lnϕV[i] + ŷ[i]) for i in 1:NC))
@NLexpression(model,Δ, ρ * (L * sL + V * sV))
@NLobjective(model,Min,ζ + Δ)

@complements(model, 0 <= V, sV >= 0) # V*sV==0
@complements(model, 0 <= L, sL >= 0) # L*sL==0


# @constraint(model, cc_sV, V * sV == 0)
# @constraint(model, cc_sV, V ⟂ sV)
# @constraint(model, cc_sL,L * sL == 0)
# @constraint(model, cc_sL, L ⟂ sL)
# @variable(model, IL4)
# @NLconstraint(model, cc_sL, AV * sL <= 1e-6)
# @NLconstraint(model, cc_sV, AL * sV <= 1e-6)

# @variable(model, keq[i in 1:NC])
# @NLconstraint(model ,def_keq[i in 1:NC], keq[i] == exp(ϕL[i] - ϕV[i]))

# @NLexpression(model, def_y[i in 1:NC], y[i] - keq[i] * x[i])
# @NLconstraint(model, def_y_con[i in 1:NC], y[i] == keq[i] * x[i])
# @NLobjective(model, Min, sum(def_y[i]^2 for i in 1:NC))
