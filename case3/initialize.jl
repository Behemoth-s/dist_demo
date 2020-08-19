function saturationT!(F, temp, pres, Vp)
    F .= Vp[:, 1] + Vp[:, 2] ./ (temp + Vp[:, 3]) + Vp[:, 4] .* log.(temp) + Vp[:, 5] .* temp.^Vp[:, 6] .- log(pres / 1000);
end

function initialize_temp(pres, Vp)
    f!(F, temp) = saturationT!(F, temp, pres, Vp)
    sat_temp0 = ones(nc) .* 200
    F0 = zeros(nc)
    r = nlsolve(f!, sat_temp0, autodiff=:forward)
    temp0 = r.zero' * x0
end

function initialize_K(pres, Temp, pc, tc, ω)
    K0 = pc ./ pres .* exp.(5.373 .* (1 .+ ω) .* (1 .- tc ./ Temp))
end

function forward_ϕL(x, pres, temp, phase, pc, tc, ω, k)
    R = 8.314e3
    m_factor = [0.480 1.574 0.176] # [0.48508, 1.55171, -0.15613]
    ac_factor = 0.42748
    b_factor = 0.08664
    m = m_factor[1] .+ m_factor[2] .* ω .- m_factor[3] .* ω.^2
    bi  = b_factor * R .* tc ./ pc
    ac = ac_factor * R^2 .* tc.^2 ./ pc
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