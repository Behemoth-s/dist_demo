function comp_flow_eq(L, x, V, y, F, z)
    return L * x + V * y - F * z
end
function ∇comp_flow_eq(g, L, x, V, y, F, z)
    g[1] = x
    g[2] = L
    g[3] = y
    g[4] = V
    g[5] = z
    g[6] = F
end
function Δ(L, V, sL, sV, ρ)
    return ρ * (L * sL + V * sV)
end
function ∇Δ(g, L, V, sL, sV, ρ)
    g[1] = ρ * sL
    g[2] = ρ * sV
    g[3] = ρ * L
    g[4] = ρ * V
end