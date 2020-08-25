include("SRK_EOS.jl")
function initialize(problem_type;z0,F0,pres0,temp0,x0,y0,L0,V0)
    if problem_type == "PT_FLASH"
        # z = z0
        # F = F0
        pres = pres0
        temp = temp0
        # x0 = z0
        if LC_EXPRESSION_FLAG
            mc0 =  F0 .* z0
            fix.(mc, mc0)
        else
            fix.(z, z0)
            fix(F, F0)
        end
        if isassigned(x0)
        end
        fix(Temp, temp0)
        set_upper_bound(L, F0)
        set_upper_bound(V, F0)

    end
    # keq0 = initialize_K(pres, temp0, PRES_CRITICAL, TEMP_CRITICAL, ω)
    # y0 = x0 .* keq0
    # y0 = y0 ./ sum(y0)


    set_start_value.(x, x0)
    set_start_value.(y, y0)



    ar0 = sqrt.(AC_VAL) .* (1 .+ MI_VAL .* (1 .- sqrt.(temp0 ./ TEMP_CRITICAL)))
    set_start_value.(ar, ar0)

    aij0 = (1 .- K) .* (ar0 * ar0')
    for i in 1:NC, j in i:NC
        set_start_value(aij[i,j], aij0[i,j])
    end

    amL0 = sum(x0[i] * x0[j] * aij0[i,j] for i in 1:NC, j in 1:NC)
    amV0 = sum(y0[i] * y0[j] * aij0[i,j] for i in 1:NC, j in 1:NC)
    set_start_value(amL, amL0)
    set_start_value(amV, amV0)

    bmL0 = sum(x0[i] * BI_VAL[i] for i in 1:NC)
    bmV0 = sum(y0[i] * BI_VAL[i] for i in 1:NC)
    set_start_value(bmL, bmL0)
    set_start_value(bmV, bmV0)

    δL0 = 2 .* ar0 ./ amL0 .* sum(x0[j] * ar0[j] .* (1 .- K[:,j]) for j in 1:NC)
    δV0 = 2 .* ar0 ./ amV0 .* sum(y0[j] * ar0[j] .* (1 .- K[:,j]) for j in 1:NC)
    set_start_value.(δL, δL0)
    set_start_value.(δV, δV0)

    AL0 = amL0 * pres / R_VAL^2 / temp0^2
    AV0 = amV0 * pres / R_VAL^2 / temp0^2
    set_start_value(AL, AL0)
    set_start_value(AV, AV0)

    BL0 = bmL0 * pres / R_VAL / temp0
    BV0 = bmV0 * pres / R_VAL / temp0
    set_start_value(BL, BL0)
    set_start_value(BV, BV0)

    ZV0 = solve_cube_Z(AV0, BV0, "V")
    ZL0 = solve_cube_Z(AL0, BL0, "L")
    set_start_value(ZL, ZL0)
    set_start_value(ZV, ZV0)

    IL10 = ZL0 - BL0;
    IV10 = ZV0 - BV0;
    set_start_value(IL1, IL10)
    set_start_value(IV1, IV10)

    IL30 = (ZL0 + BL0) / ZL0;
    IV30 = (ZV0 + BV0) / ZV0;
    set_start_value(IL3, IL30)
    set_start_value(IV3, IV30)

    if LN_EXPRESSION_PHI
        set_start_value.(x̂, log.(x0))
        set_start_value.(ŷ, log.(y0))
        lnϕL0 = (ZL0 - 1) .* BI_VAL ./ bmL0 .- IL10 .- AL0 ./ BL0 .* (δL0  .- BI_VAL ./ bmL0) .* IL30;
        lnϕV0 = (ZV0 - 1) .* BI_VAL ./ bmV0 .- IV10 .- AV0 ./ BV0 .* (δV0  .- BI_VAL ./ bmV0) .* IV30;
        set_start_value.(lnϕL, lnϕL0)
        set_start_value.(lnϕV, lnϕV0)
        det0 = sum(y0 .- exp.(lnϕL0 .- lnϕV0) .* x0)
        set_lower_bound.(x̂, ϵx)
        set_upper_bound.(x̂, 0)
        set_lower_bound.(ŷ, ϵx)
        set_upper_bound.(ŷ, 0)
        set_start_value(β, 0)


    else
        ϕL0 = exp.((ZL0 - 1) .* BI_VAL ./ bmL0 .- IL10 .- AL0 ./ BL0 .* (δL0  .- BI_VAL ./ bmL0) .* IL30);
        ϕV0 = exp.((ZV0 - 1) .* BI_VAL ./ bmV0 .- IV10 .- AV0 ./ BV0 .* (δV0  .- BI_VAL ./ bmV0) .* IV30);
        set_start_value.(ϕL, ϕL0)
        set_start_value.(ϕV, ϕV0)
        det0 = sum(y0 .- ϕL0 ./ ϕL0 .* x0)
        set_start_value(β, 0)

    end

    set_start_value(sL, 0.000001)
    set_start_value(sV, 0.000001)

    set_start_value(L, L0)
    set_start_value(V, V0)

    if LC_EXPRESSION_FLAG
        set_start_value.(lc, L0 .* x0)
        set_start_value.(vc, V0 .* y0)
    end
    set_lower_bound.(x, 0)
    set_upper_bound.(x, 1)
    set_lower_bound.(y, 0)
    set_upper_bound.(y, 1)
    det1 = sum(F0 .* z0 - L0 .* x0 - V0 .* y0)

    println("y-kx", det0)
    println("F-LV", det1)



end

function initlize_bubleT(x0, pres)

    temp0 = initialize_temp(pres, VAP_PRE)
    keq0 = initialize_K(pres, temp0, PRES_CRITICAL, TEMP_CRITICAL, ω)
    y0 = x0 .* keq0
    y0 = y0 ./ sum(y0)



    ar0 = sqrt.(AC_VAL) .* (1 .+ MI_VAL .* (1 .- sqrt.(temp0 ./ TEMP_CRITICAL)))
    aij0 = (1 .- K) .* (ar0 * ar0')
    amL0 = sum(x0[i] * x0[j] * aij0[i,j] for i in 1:NC, j in 1:NC)
    bmL0 = sum(x0[i] * BI_VAL[i] for i in 1:NC)
    AL0 = amL0 * pres / R_VAL^2 / temp0^2
    BL0 = bmL0 * pres / R_VAL / temp0
    ZL0 = solve_cube_Z(AL0, BL0, "L")
    atmpL0 = (1 .- K) * (x0 .* ar0);
    lnϕL0 = (ZL0 - 1) .* BI_VAL / bmL0 .- log(ZL0 - BL0) .- AL0 ./ BL0 .* (2 .* ar0 ./ amL0 .* atmpL0 .- BI_VAL ./ bmL0) .* log(1 + BL0 / ZL0);

    amV0 = sum(y0[i] * y0[j] * aij0[i,j] for i in 1:NC, j in 1:NC)
    bmV0 = sum(y0[i] * BI_VAL[i] for i in 1:NC)
    AV0 = amV0 * pres / R_VAL^2 / temp0^2
    BV0 = bmV0 * pres / R_VAL / temp0
    ZV0 = solve_cube_Z(AV0, BV0, "V")
    atmpV0 = (1 .- K) * (y0 .* ar0);
    lnϕV0 = (ZV0 - 1) .* BI_VAL / bmV0 .- log(ZV0 - BV0) .- AV0 ./ BV0 .* (2 .* ar0 ./ amV0 .* atmpV0 .- BI_VAL ./ bmV0) .* log(1 + BV0 / ZV0);


    fix.(x, x0)
    set_start_value.(y, y0)
    # set_start_value.(keq, keq0)
    set_start_value(Temp, temp0)
    set_start_value(AL, AL0)
    set_start_value(BL, BL0)
    set_start_value(ZL, ZL0)
    set_start_value(AV, AV0)
    set_start_value(BV, BV0)
    set_start_value(ZV, ZV0)

    set_start_value.(ar, ar0)
    for i in 1:NC, j in i:NC
        set_start_value(aij[i, j], aij0[i,j])
    end
    set_start_value(amL, amL0)
    set_start_value(amV, amV0)

    set_start_value.(atmpL, atmpL0)
    set_start_value.(atmpV, atmpV0)

    set_start_value(bmL, bmL0)
    set_start_value(bmV, bmV0)
    set_start_value.(ϕL, exp.(lnϕL0))
    set_start_value.(ϕV, exp.(lnϕV0))
    set_start_value(sL, 0)
    set_start_value(sV, 0)
    set_start_value(β, 1)
end
function saturationT!(F, temp, pres, Vp)
    F .= Vp[:, 1] + Vp[:, 2] ./ (temp + Vp[:, 3]) + Vp[:, 4] .* log.(temp) + Vp[:, 5] .* temp.^Vp[:, 6] .- log(pres / 1000);
end

function initialize_temp(pres, Vp)
    nc  = size(Vp)[1]
    f!(F, temp) = saturationT!(F, temp, pres, Vp)
    sat_temp0 = ones(nc) .* 200
    F0 = zeros(nc)
    r = nlsolve(f!, sat_temp0, autodiff=:forward)
    temp0 = r.zero' * x0
end

function initialize_K(pres, temp, pc, tc, ω)
    K0 = pc ./ pres .* exp.(5.373 .* (1 .+ ω) .* (1 .- tc ./ temp))
end




