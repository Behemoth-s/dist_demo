using JuMP
using MathOptInterface
using Complementarity
const R_VAL = 8.314e3
const M_FACTOR = [0.480 1.574 0.176] # [0.48508, 1.55171, -0.15613]
const AC_FACTOR = 0.42748
const B_FACTOR = 0.08664
const MI_VAL = M_FACTOR[1] .+ M_FACTOR[2] .* ω .- M_FACTOR[3] .* ω.^2
const BI_VAL  = B_FACTOR * R_VAL .* TEMP_CRITICAL ./ PRES_CRITICAL
const AC_VAL = AC_FACTOR * R_VAL^2 .* TEMP_CRITICAL.^2 ./ PRES_CRITICAL
const u = 1
const v = 0
const I0 = sqrt(u^2 - 4v)
const K_MEAN = sum(K) / NC / NC
if USER_FUNCTION_FLAG
    include("SRK_EOS.jl")
    register(model, :cube_eos_Z, 3, SRK_EOS_Z, ∇SRK_EOS_Z)
    register(model, :cube_eos_dZ, 3, SRK_EOS_dZ, ∇SRK_EOS_dZ)
    register(model, :cube_eos_lnϕ, 9, SRK_EOS_lnϕ, ∇SRK_EOS_lnϕ)
    register(model, :cube_eos_ϕ, 9, SRK_EOS_ϕ, ∇SRK_EOS_ϕ)
    register(model, :cube_eos_A, 4, SRK_EOS_A, ∇SRK_EOS_A)
    register(model, :cube_eos_B, 4, SRK_EOS_B, ∇SRK_EOS_B)
    register(model, :cube_eos_I1, 3, SRK_EOS_I1, ∇SRK_EOS_I1)
    register(model, :cube_eos_I3, 3, SRK_EOS_I3, ∇SRK_EOS_I3)
    include("FLASH_EO.jl")
    register(model, :comp_flow_eq, 6, comp_flow_eq, ∇comp_flow_eq)
    register(model, :Δ, 5, Δ, ∇Δ)
end
@NLparameter(model, bi_val[i in 1:NC] == BI_VAL[i])
@variable(model, Temp >= 0)
@variable(model, Pres >= 0)

@variable(model, x[1:NC])
@variable(model, y[1:NC])

@variable(model, IL1)
@variable(model, IV1)

@variable(model, IL3)
@variable(model, IV3)

@variable(model, ar[1:NC] >= 0)
@variable(model, aij[i in 1:NC, j in i:NC] >= 0)

@variable(model, amL >= 0)
@variable(model, amV >= 0)

@variable(model, δL[1:NC])
@variable(model, δV[1:NC])

@variable(model, bmL )
@variable(model, bmV )

@variable(model, AV )
@variable(model, AL )

@variable(model, BV )
@variable(model, BL )

@variable(model, ZV >= 0)
@variable(model, ZL >= 0)




@variable(model, β)

@NLconstraint(model, def_ar[i in 1:NC], ar[i] == sqrt(AC_VAL[i]) * (1 + MI_VAL[i] * (1 - sqrt(Temp / TEMP_CRITICAL[i]))))
@constraint(model, def_aij[i in 1:NC, j in i:NC], aij[i, j] == (1 - K[i,j]) * ar[i] * ar[j])

@NLconstraint(model, def_amL, amL == sum(x[i]^2 * aij[i,i] for i in 1:NC) + 2 * sum(x[i] * x[j] * aij[i,j] for i in 1:NC, j in i + 1:NC))
@NLconstraint(model, def_amV, amV == sum(y[i]^2 * aij[i,i] for i in 1:NC) + 2 * sum(y[i]^2 * aij[i,j] for i in 1:NC, j in i + 1:NC))
if K_MEAN_FLAG
    @NLconstraint(model, def_δL[i in 1:NC], δL[i] == 2 * ar[i] / amL * sum(x[j] * ar[j] * (1 - K_MEAN) for j in 1:NC))
    @NLconstraint(model, def_δV[i in 1:NC], δV[i] == 2 * ar[i] / amV * sum(y[j] * ar[j] * (1 - K_MEAN) for j in 1:NC))
else
    @NLconstraint(model, def_δL[i in 1:NC], δL[i] == 2 * ar[i] / amL * sum(x[j] * ar[j] * (1 - K[i,j]) for j in 1:NC))
    @NLconstraint(model, def_δV[i in 1:NC], δV[i] == 2 * ar[i] / amV * sum(y[j] * ar[j] * (1 - K[i,j]) for j in 1:NC))
end
@constraint(model, def_bmV, bmV == sum(y[i] * BI_VAL[i] for i in 1:NC))
@constraint(model, def_bmL, bmL == sum(x[i] * BI_VAL[i] for i in 1:NC))

# if LN_EXPRESSION_A
#     @NLconstraint(model, def_AL, log(AL) == log(amL) + log(Pres) - 2 * log(R_VAL) - 2log(Temp))
#     @NLconstraint(model, def_AV, log(AV) == log(amV) + log(Pres) - 2 * log(R_VAL) - 2log(Temp))
#     @NLconstraint(model, def_BV, log(BV) == log(bmV) + log(Pres) - log(R_VAL) - log(Temp))
#     @NLconstraint(model, def_BL, log(BL) == log(bmL) + log(Pres) - log(R_VAL) + log(Temp))

# else
if USER_FUNCTION_FLAG
    @NLconstraint(model, def_AL, cube_eos_A(AL, amL, Temp, Pres) == 0)
    @NLconstraint(model, def_AV, cube_eos_A(AV, amV, Temp, Pres) == 0)
    @NLconstraint(model, def_BL, cube_eos_B(BL, bmL, Temp, Pres) == 0)
    @NLconstraint(model, def_BV, cube_eos_B(BV, bmV, Temp, Pres) == 0)

else
    @NLconstraint(model, def_AL, AL == amL * Pres / R_VAL^2 / Temp^2)
    @NLconstraint(model, def_BL, BL == bmL * Pres / R_VAL / Temp)
    @NLconstraint(model, def_AV, AV == amV * Pres / R_VAL^2 / Temp^2)
    @NLconstraint(model, def_BV, BV == bmV * Pres / R_VAL / Temp)
end
# end

# 简化lnϕ
if USER_FUNCTION_FLAG
    @NLconstraint(model, def_IL1,cube_eos_I1(IL1, ZL, BL) == 0)
    @NLconstraint(model, def_IV1,cube_eos_I1(IV1, ZV, BV) == 0)
    @NLconstraint(model, def_IL3,cube_eos_I3(IL3, ZL, BL) == 0)
    @NLconstraint(model, def_IV3,cube_eos_I3(IV3, ZV, BV) == 0)
else
    @NLconstraint(model,def_IL1, exp(IL1) == ZL - BL)
    @NLconstraint(model,def_IV1, exp(IV1) == ZV - BV)

    @NLconstraint(model, def_IL3, exp(IL3) == ( ZL + BL ) / ZL)
    @NLconstraint(model, def_IV3, exp(IV3) == ( ZV + BV ) / ZV)
end


if LN_EXPRESSION_PHI
    @variable(model, lnϕV[1:NC])
    @variable(model, lnϕL[1:NC])
    @variable(model, x̂[1:NC])
    @variable(model, ŷ[1:NC])
    if USER_FUNCTION_FLAG
        @NLconstraint(model, def_phiL[i in 1:NC], cube_eos_lnϕ(lnϕL[i], ZL, AL, BL, bmL, δL[i], IL1, IL3, bi_val[i]) == 0)
        @NLconstraint(model, def_phiV[i in 1:NC], cube_eos_lnϕ(lnϕV[i], ZV, AV, BV, bmV, δV[i], IV1, IV3, bi_val[i]) == 0)
    else
        @NLconstraint(model, def_phiL[i in 1:NC], lnϕL[i] == (ZL - 1) * BI_VAL[i] / bmL - IL1 - AL / BL * (δL[i]  - bi_val[i] / bmL) * IL3)
        @NLconstraint(model, def_phiV[i in 1:NC], lnϕV[i] == (ZV - 1) * BI_VAL[i] / bmV - IV1 - AV / BV * (δV[i]  - bi_val[i] / bmV) * IV3)
    end
    @NLconstraint(model, def_lnx[i in 1:NC], x[i] == exp(x̂[i]))
    @NLconstraint(model, def_lny[i in 1:NC], y[i] == exp(ŷ[i]))
    @NLconstraint(model, phase_eq[i in 1:NC], ŷ[i] + lnϕV[i] == β + lnϕL[i] + x̂[i])

else
    @variable(model, ϕV[1:NC])
    @variable(model, ϕL[1:NC])
    if USER_FUNCTION_FLAG
        @NLconstraint(model, def_phiL[i in 1:NC], cube_eos_ϕ(ϕL[i], ZL, AL, BL, bmL, δL[i], IL1, IL3, bi_val[i]) == 0)
        @NLconstraint(model, def_phiV[i in 1:NC], cube_eos_ϕ(ϕV[i], ZV, AV, BV, bmV, δV[i], IV1, IV3, bi_val[i]) == 0)
    else
        @NLconstraint(model, def_phiL[i in 1:NC], log(ϕL[i]) == ((ZL - 1) * BI_VAL[i] / bmL - IL1 - AL / BL * (δL[i]  - bi_val[i] / bmL) * IL3))
        @NLconstraint(model, def_phiV[i in 1:NC], log(ϕV[i]) == ((ZV - 1) * BI_VAL[i] / bmV - IV1 - AV / BV * (δV[i]  - bi_val[i] / bmV) * IV3))
    end
    @NLconstraint(model, phase_eq[i in 1:NC], y[i] == β * ϕL[i] / ϕV[i] * x[i] )

end
# 确定三次方程根
if USER_FUNCTION_FLAG
    @NLconstraint(model, def_ZL1, cube_eos_Z(ZL, AL, BL) == 0)
    @NLconstraint(model, def_ZV1, cube_eos_Z(ZV, AV, BV) == 0)
    # @NLconstraint(model, def_ZL2, cube_eos_dZ(ZL, AL, BL) >= 0)
    # @NLconstraint(model, def_ZV2, cube_eos_dZ(ZV, AV, BV) >= 0)

else
    @NLconstraint(model, def_ZL1, ZL^3 - ZL^2 + (AL - BL^2 - BL) * ZL - AL * BL == 0)
    @NLconstraint(model, def_ZV1, ZV^3 - ZV^2 + (AV - BV^2 - BV) * ZV - AV * BV == 0)
    @NLconstraint(model, def_ZL2, 3 * ZL^2 - 2 *  ZL + (AL - BL^2 - BL) >= 0)
    @NLconstraint(model, def_ZV2, 3 * ZV^2 - 2 * ZV + (AV - BV^2 - BV) >= 0)
end




if PROBLEM_TYPE == "PT_FLASH"
    @variable(model, L >= 0)
    @variable(model, V >= 0)
    @variable(model,sL >= 0)
    @variable(model,sV >= 0)
    if LN_EXPRESSION_PHI
        @constraint(model, lower_limit_β, β  >= -sL)
        @constraint(model, upper_limit_β, β  <= sV)
    else
        @constraint(model, lower_limit_β, β - 1  >= -sL)
        @constraint(model, upper_limit_β, β - 1  <= sV)
    end
    if COMPLEMENTS_FLAG
        @complements(model, 0 <= V, sV >= 0) # V*sV==0
        @complements(model, 0 <= L, sL >= 0) # L*sL==0
    else
        @constraint(model, V * sV == 0)
        @constraint(model, L * sL == 0)
    end

    if LC_EXPRESSION_FLAG
        @variable(model, lc[1:NC])
        @variable(model, vc[1:NC])
        @variable(model, mc[1:NC])
        @constraint(model, lv_eq[i in 1:NC], mc[i] == lc[i] + vc[i])

        @constraint(model, def_l[i in 1:NC], lc[i] == L * x[i])
        @constraint(model, def_v[i in 1:NC], vc[i] == V * y[i])

        @constraint(model, sum_l, L == sum(lc[i] for i in 1:NC))
        @constraint(model, sum_v, V == sum(vc[i] for i in 1:NC))
    else
        @variable(model, z[1:NC])
        @variable(model, F)
        @constraint(model, eq_flow, F == L + V)
        if USER_FUNCTION_FLAG
            @constraint(model, eq_cmp[i in 1:NC], comp_flow_eq(L, x[i], V, y[i], F, z[i]) == 0)
        else
            @constraint(model, eq_cmp[i in 1:NC], F * z[i] == L * x[i] + V * y[i])
        end
    end
    @constraint(model, sumy, sum(y[i] - x[i] for i in 1:NC) == 0 )
    @NLobjective(model,Min, sum((ŷ[i] + lnϕV[i] -  β - lnϕL[i] - x̂[i])^2 for i in 1:NC))
    @constraint(model, def_ZL3, 6 * ZL - 2   <= M_COFF * sL)
    @constraint(model, def_ZV3, 6 * ZV - 2   >= -M_COFF * sV)

elseif PROBLEM_TYPE == "BUBBLET"

    # @variable(model,sL >= 0)
    # @variable(model,sV >= 0)
    # set_start_value(sL, 0)
    # set_start_value(sV, 0)
    # @constraint(model, def_ZL3, 6 * ZL - 2   <= 0)
    # @constraint(model, def_ZV3, 6 * ZV - 2   >= 0)

    # if LN_EXPRESSION_PHI
    #     @constraint(model, lower_limit_β, β  >= -sL)
    #     @constraint(model, upper_limit_β, β  <= sL)
    # else
    #     @constraint(model, lower_limit_β, β - 1  >= -sL)
    #     @constraint(model, upper_limit_β, β - 1  <= sL)
    # end
    if LN_EXPRESSION_PHI
        @constraint(model, complements_temp, Temp * β == 0)
        @NLobjective(model, Min, 1e6 * (sum(y[i] for i in 1:NC) - 1)^2 + Temp * (β))
    else
        # @constraint(model, complements_temp, Temp * (β - 1) == 0)

        # @variable(model, sy >= 0)
        # @variable(model, sx >= 0)

        # @constraint(model, defsumy, sum(y[i] for i in 1:NC) - 1 == 0)
        # @constraint(model, defsumx, sum(y[i] for i in 1:NC) - 1 >= -sx)
        @NLconstraint(model, sumy, sum(ϕL[i] / ϕV[i] * x[i] for i in 1:NC) - 1 == 0)
        # @variable(model,sL >= 0)
        # @constraint(model, lower_limit_β, β - 1  >= -sy)
        # @constraint(model, upper_limit_β, β - 1  <= sx)
        # @NLobjective(model, Min, (sum(y[i] for i in 1:NC) - 1)^2)
        @variable(model, sP >= 0)
        @constraint(model, sP >= PRES - Pres)
        @constraint(model, sP >= Pres - PRES )
        set_lower_bound.(y, 0)
        set_upper_bound.(y, 1)
        # set_lower_bound(β, 0.7)
        set_start_value(Pres, PRES - 1)
        @NLobjective(model, Min, 1e-3 * sP * Temp)
        fix(β, 1)
        # @NLobjective(model, Min,  1e6 * (sum(y[i] for i in 1:NC) - 1)^2)
    end
end









# @NLexpression(model, ζ, sum(lc[i] * (lnϕL[i] + x̂[i]) + vc[i] * (lnϕV[i] + ŷ[i]) for i in 1:NC))

# @NLexpression(model,Δ, ρ * (L * sL + V * sV))



