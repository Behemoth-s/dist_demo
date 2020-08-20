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
    set_lower_bound(ZL, 0)
    set_lower_bound(ZV, 0)
    set_start_value(sL, 0)
    set_start_value(sV, 0)
    set_start_value(β, 1)
    set_lower_bound.(y, 0)
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


function forward_ϕL(x, pres, temp, phase, pc, tc, ω, k)

    m = M_FACTOR[1] .+ M_FACTOR[2] .* ω .- M_FACTOR[3] .* ω.^2
    bi  = B_FACTOR * R .* tc ./ pc
    ac = AC_FACTOR * R^2 .* tc.^2 ./ pc
    ar = sqrt.(ac) .* (1 .+ m .* (1 .- sqrt.(temp ./ tc)))
    aij = (1 .- k) .* (ar * ar')
    am = sum(x[i] * x[j] * aij[i,j] for i in 1:nc, j in 1:nc)
    bm = sum(x[i] * bi[i] for i in 1:nc)
    A = am * pres / R^2 / temp^2
    B = bm * pres / R / temp
    Z = solve_cube_Z(A, B, phase)
    atmp = (1 .- k) * (x .* ar);
    lnphi = (Z - 1) .* bi / bm .- log(abs(Z - B)) .- A ./ B .* (2 .* ar ./ am .* atmp .- bi ./ bm) .* log(1 + B / Z);
    # ϕ = (Z - 1) .* bi ./ bm .- log(Z - B) .- A ./ B .* (2 .* ami ./ am  .- bi ./ bm) .* log(1 + B / Z)
end

function solve_cube_Z(A, B, phase="L")
    poly = zeros(4)
    poly[4] = 1
    poly[3] = -1
    poly[2] = A - B^2 - B
    poly[1] = -A * B
    roots, flag = solve_cubic_eq(float.(complex(poly)))
    roots_real = real.(roots)
    roots_imag = imag.(roots)
    Z = 1
    if flag == 0
        Z = real(roots[1])
    end
    if flag < 0
        if phase == "L"
            Z = min(real_roots[real_roots .> 0])
        else
            Z = max(real_roots[real_roots .> 0])
        end

    end
    if flag > 0
        if roots_imag[1] + roots_imag[2] ≈ 0
            Z =  roots_real[3]
        elseif roots_imag[1] + roots_imag[3] ≈ 0
            Z =  roots_real[2]
        else
            Z =  roots_real[1]
        end
    end
    return Z

end

function solve_cubic_eq(poly::AbstractVector{Complex{T}}) where {T <: AbstractFloat}
# Cubic equation solver for complex polynomial (degree=3)
# http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
    third = 1 // 3
    a1  =  1 / poly[4]
    E1  = -poly[3] * a1
    E2  =  poly[2] * a1
    E3  = -poly[1] * a1
    s0  =  E1
    E12 =  E1 * E1
    A   =  2 * E1 * E12 - 9 * E1 * E2 + 27 * E3 # = s1^3 + s2^3
    B   =  E12 - 3 * E2                 # = s1 s2
    # quadratic equation: z^2 - Az + B^3=0  where roots are equal to s1^3 and s2^3
    Δ = sqrt(A * A - 4 * B * B * B)
    if real(conj(A) * Δ) >= 0 # scalar product to decide the sign yielding bigger magnitude
        s1 = exp(log(0.5 * (A + Δ)) * third)
    else
        s1 = exp(log(0.5 * (A - Δ)) * third)
    end
    if s1 == 0
        s2 = s1
    else
        s2 = B / s1
    end
    zeta1 = complex(-0.5, sqrt(T(3.0)) * 0.5)
    zeta2 = conj(zeta1)
    return (third * (s0 + s1 + s2), third * (s0 + s1 * zeta2 + s2 * zeta1), third * (s0 + s1 * zeta1 + s2 * zeta2)), sign(real(Δ))
end