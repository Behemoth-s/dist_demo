function saturationT!(F, temp, pres, Vp)
    F .= Vp[:, 1] + Vp[:, 2] ./ (temp + Vp[:, 3]) + Vp[:, 4] .* log.(temp) + Vp[:, 5] .* temp.^Vp[:, 6] .- log(pres / 1000);
end
function initlize_temp(pres, Vp)
    f!(F, temp) = saturationT!(F, temp, pres, Vp)
    sat_temp0 = ones(nc) .* 200
    F0 = zeros(nc)
    r = nlsolve(f!, sat_temp0, autodiff=:forward)
    temp0 = r.zero' * x0
end
function initlize_K(pres, Temp, pc, Tc, ω)
    K0 = Pc ./ pres .* exp.(5.373 .* (1 .+ ω) .* (1 .- Tc ./ Temp))
end