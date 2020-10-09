
function SRK_EOS_a(a, Temp, ac_val, mi_val, tc)
    return a -  ac_val * (1 + mi_val * (1 - sqrt(Temp / tc)))
end
function ∇SRK_EOS_a(g, a, Temp, ac_val, mi_val, tc)
    g[1] = 1
    g[2] = -ac_val * mi_val^2 / tc +  ac_val * (1 + mi_val) * sqrt(tc / Temp)
end
function SRK_EOS_aij(aij, ai, aj, factor_k)
    """math
    a_{ij} = (1 - K_{ij})
    """
    # a_ij = (1 - K[i,j]) * ar[i] * ar[j]
    return log((1 - factor_k)) + log(ai) / 2 + log(aj) / 2 - log(aij)
end
function ∇SRK_EOS_aij(g, aij, ai, aj, factor_k)
    g[1] = -1 / aij
    g[2] = 1 / 2 / ai
    g[3] = 1 / 2 / aj
end



function SRK_EOS_lnϕ(lnϕ, Z, A, B, bm, δ, I1, I3, bi)
    return lnϕ -  (Z - 1) * bi / bm + I1 + A / B * (δ  - bi / bm) * I3
end
function ∇SRK_EOS_lnϕ(g, lnϕ, Z, A, B, bi, bm, δ, I1, I3)
    # lnϕ
    g[1] = 1
    # Z
    g[2] = - bi / bm
    # A
    g[3] = (δ - bi / bm) * I3 / B
    # B
    g[4] = -A * (δ  - bi / bm) * I3 / B^2
    # bm
    g[5] = ((Z - 1)  + A / B * I3) * bi / bm^2
    # δ
    g[6] = A / B * I3
    # I1
    g[7] = 1
    # I3
    g[8] = A / B * (δ  - bi / bm)
# bi:parameter
    g[9] = -(Z - 1) - A / B * I3 / bm
end

function SRK_EOS_Z(Z, A, B)
    return Z^3 - Z^2 + (A - B^2 - B) * Z - A * B
end
function ∇SRK_EOS_Z(g, Z, A, B)
    g[1] = 3 * Z^2 - 2 * Z + (A - B^2 - B)
    g[2] = Z - B
    g[3] = -2 * Z * B - Z - A
end

function SRK_EOS_dZ(Z, A, B)
    return 3 * Z^2 - 2 *  Z + (A - B^2 - B)
end
function ∇SRK_EOS_dZ(g, Z, A, B)
    g[1] = 6 * Z - 2
    g[2] = 1
    g[3] = -2 * B - 1
end

function SRK_EOS_A(A, am, Temp)
    return A - am * PRES / R_VAL^2 / Temp^2
end
function ∇SRK_EOS_A(g, A, am, Temp)
    g[1] = 1
    g[2] = -PRES / R_VAL^2 / Temp^2
    g[3] = 2 * am * PRES / R_VAL^2 / Temp^3
end

function SRK_EOS_B(B, bm, Temp)
    return B - bm * PRES / R_VAL / Temp
end
function ∇SRK_EOS_B(g, B, bm, Temp)
    g[1] = 1
    g[2] = -PRES / R_VAL / Temp
    g[3] = bm * PRES / R_VAL / Temp^2
end

function SRK_EOS_I1(I1, Z, B)
    return exp(I1) - Z + B
end
function ∇SRK_EOS_I1(g, I1, Z, B)
    g[1] = exp(I1)
    g[2] = -1
    g[3] = 1
end
function SRK_EOS_I3(I3, Z, B)
    return exp(I3) - 1 - B / Z
end
function ∇SRK_EOS_I3(g, I3, Z, B)
    g[1] = exp(I3)
    g[2] = B / Z^2
    g[3] = -1 / Z
end

function forward_lnϕ(z::Array{Float64,1}, pres::Float64, temp::Float64, liq_phase::Bool)::Float64
    # cacluate lnϕ
    # m = M_FACTOR[1] .+ M_FACTOR[2] .* ω .- M_FACTOR[3] .* ω.^2
    # bi  = B_FACTOR * R .* tc ./ pc
    # ac = AC_FACTOR * R^2 .* tc.^2 ./ pc
    ar = sqrt.(AC_VAL) .* (1 .+ MI_VAL .* (1 .- sqrt.(temp ./ TEMP_CRITICAL)))
    aij = (1 .- K) .* (ar * ar')
    am = sum(z[i] * z[j] * aij[i,j] for i in 1:nc, j in 1:nc)
    bm = sum(z[i] * BI_VAL[i] for i in 1:nc)
    A = am * pres / R^2 / temp^2
    B = bm * pres / R / temp
    Z = solve_cube_Z(A, B, liq_phase)
    atmp = (1 .- k) * (z .* ar)
    return (Z - 1) .* bi / bm .- log(abs(Z - B)) .- A ./ B .* (2 .* ar ./ am .* atmp .- bi ./ bm) .* log(1 + B / Z)
    # ϕ = (Z - 1) .* bi ./ bm .- log(Z - B) .- A ./ B .* (2 .* ami ./ am  .- bi ./ bm) .* log(1 + B / Z)
end

function solve_cube_Z(A::Float64, B::Float64, liq_phase::Bool=true)::Float64
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
        if liq_phase
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