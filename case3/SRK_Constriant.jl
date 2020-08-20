using JuMP
using MathOptInterface
const R_VAL = 8.314e3
M_FACTOR = [0.480 1.574 0.176] # [0.48508, 1.55171, -0.15613]
const AC_FACTOR = 0.42747
const B_FACTOR = 0.08664
const MI_VAL = M_FACTOR[1] .+ M_FACTOR[2] .* ω .- M_FACTOR[3] .* ω.^2
const BI_VAL  = B_FACTOR * R_VAL .* TEMP_CRITICAL ./ PRES_CRITICAL
const AC_VAL = AC_FACTOR * R_VAL^2 .* TEMP_CRITICAL.^2 ./ PRES_CRITICAL
const M_COFF = 10

@variable(model, AL)
@variable(model, BL)
@variable(model, ZL)

@variable(model, ar[1:NC])
@variable(model, aij[i in 1:NC, j in i:NC])

@variable(model, amL)
@variable(model, atmpL[1:NC])

@variable(model, bmL)
@variable(model, ϕL[i in 1:NC])
@variable(model, AV)
@variable(model, BV)
@variable(model, ZV)

@variable(model, amV)
@variable(model, atmpV[1:NC])

@variable(model, bmV)
@variable(model, ϕV[i in 1:NC])
@variable(model, Temp)
@variable(model, x[i in 1:NC])
@variable(model, y[1:NC])

@NLconstraint(model, def_ar[i in 1:NC], ar[i] == sqrt(AC_VAL[i]) * (1 + MI_VAL[i] * (1 - sqrt(Temp / TEMP_CRITICAL[i]))))
@NLconstraint(model, def_aij[i in 1:NC, j in i:NC], aij[i, j] == (1 - K[i,j]) * ar[i] * ar[j])

@NLconstraint(model, def_amL, amL == sum(x[i]^2 * aij[i,i] for i in 1:NC) + 2 * sum(x[i] * x[j] * aij[i,j] for i in 1:NC, j in i + 1:NC))
@NLconstraint(model, def_atmpL[i in 1:NC], atmpL[i] == sum((1 - K[i, j]) * (x[j] * ar[j]) for j in 1:NC))
# @NLconstraint(model, def_amiL[i in 1:nc], amiL[i] / x[i] == sum(aij[j, i] for j in 1:i - 1) + sum(aij[i, j] for j in i:nc))
@NLconstraint(model, def_bmL, bmL == sum(x[i] * BI_VAL[i] for i in 1:NC))

@NLconstraint(model, def_AL, AL == amL * PRES / R_VAL^2 / Temp^2)
@NLconstraint(model, def_BL, BL == bmL * PRES / R_VAL / Temp)

@NLconstraint(model, def_phiL[i in 1:NC], ϕL[i] == exp((ZL - 1) * BI_VAL[i] / bmL - log(ZL - BL) - AL / BL * (2 * ar[i] / amL * atmpL[i] - BI_VAL[i] / bmL) * log(1 + BL / ZL)))

# @NLconstraint(model, def_phiL[i in 1:nc], ϕL[i] == (ZL - 1) * bi[i] / bmL - log(ZL - BL) - AL / BL * (2 * amiL[i] / amL  - bi[i] / bmL) * log(1 + BL / ZL))

@NLconstraint(model, def_ZL1, ZL^3 - ZL^2 + (AL - BL^2 - BL) * ZL - AL * BL == 0)
@NLconstraint(model, def_ZL2, 3 * ZL^2 - 2 *  ZL + (AL - BL^2 - BL) >= 0)
# @NLconstraint(model, def_ZL3, 6 * ZL - 2   <= 0)







# @NLconstraint(model, def_arV[i in 1:nc], arV[i] == sqrt(ac[i]) * (1 + m[i] * (1 - Temp / Tc[i])))
# @NLconstraint(model, def_aij[i in 1:nc, j in i:nc], aij[i, j] == (1 - k[i,j]) * ar[i] * ar[j])
@NLconstraint(model, def_amV, amV == sum(y[i]^2 * aij[i,i] for i in 1:NC) + 2 * sum(y[i]^2 * aij[i,j] for i in 1:NC, j in i + 1:NC))
@NLconstraint(model, def_atmpV[i in 1:NC], atmpV[i] == sum((1 - K[i, j]) * (y[j] * ar[j]) for j in 1:NC))

# @NLconstraint(model, def_amiV[i in 1:nc], amiV[i] / y[i] == sum(aij[j, i] for j in 1:i - 1) + sum(aij[i, j] for j in i:nc))
@NLconstraint(model, def_bmV, bmV == sum(y[i] * BI_VAL[i] for i in 1:NC))
@NLconstraint(model, def_AV, AV == amV * PRES / R_VAL^2 / Temp^2)
@NLconstraint(model, def_BV, BV == bmV * PRES / R_VAL^2 / Temp^2)

# @NLconstraint(model, def_phiV[i in 1:nc], ϕV[i] == (ZV - 1) * bi[i] / bmV - log(ZV - BV) - AV / BV * (2 * amiV[i] / amV  - bi[i] / bmV) * log(1 + BV / ZV))
# @NLconstraint(model, def_phiV[i in 1:NC], ϕV[i] == (ZV - 1) * BI_VAL[i] / bmV - log(ZV - BV) - AV / BV * (2 * ar[i] / amV * atmpV[i] - BI_VAL[i] / bmV) * log(1 + BL / ZL))
@NLconstraint(model, def_phiV[i in 1:NC], ϕV[i] == exp((ZV - 1) * BI_VAL[i] / bmV - log(ZV - BV) - AV / BV * (2 * ar[i] / amV * atmpV[i] - BI_VAL[i] / bmV) * log(1 + BL / ZL)))

@NLconstraint(model, def_ZV1, ZV^3 - ZV^2 + (AV - BV^2 - BV) * ZV - AV * BV == 0)
@NLconstraint(model, def_ZV2, 3 * ZV^2 - 2 * ZV + (AV - BV^2 - BV) >= 0)

# @NLconstraint(model, def_ZV3, 6 * ZV - 2   >= 0)
# @NLconstraint(model, def_sumy, sum(y[i] for i in 1:NC) == 1)
@variable(model,sL >= 0)
@variable(model,sV >= 0)
@variable(model, β)
@constraint(model, lower_limit_β, β - 1 >= -sL)
@constraint(model, upper_limit_β, β - 1 <= sV)
@NLconstraint(model, def_ZV3, 6 * ZV - 2   <= -M_COFF * sV)
@NLconstraint(model, def_ZL3, 6 * ZL - 2   <= M_COFF * sL)
@NLconstraint(model, def_y[i in 1:NC], y[i] == β * ϕL[i] / ϕV[i] * x[i])
@NLconstraint(model, cc_sL, AV * sL <= 1e-6)
@NLconstraint(model, cc_sV, AL * sV <= 1e-6)

# @variable(model, keq[i in 1:NC])
# @NLconstraint(model ,def_keq[i in 1:NC], keq[i] == exp(ϕL[i] - ϕV[i]))

# @NLexpression(model, def_y[i in 1:NC], y[i] - keq[i] * x[i])
# @NLconstraint(model, def_y_con[i in 1:NC], y[i] == keq[i] * x[i])
# @NLobjective(model, Min, sum(def_y[i]^2 for i in 1:NC))
