using MacroModelling, StatsPlots
using SymPy

@model rbc begin
    γ * c[0]^-1 = γ * β * c[1]^-1 * (1 - δ + r[1])
    w[0] = - (-φ * (1 - n[0])^-1) / (γ * c[0]^-1)
    k[0] = (1-δ) * k[-1] + iv[0]
    y[0] = c[0] + iv[0]
    y[0] = a[0] * k[-1]^α * n[0]^(1-α)
    w[0] = (1-α) * y[0]/n[0]
    r[0] = α * y[0] / k[-1]
    log(a[0])= (1-ρ) * log(A) + ρ * log(a[-1]) + eps_a[x]
end

@parameters rbc begin
    β = 0.99
    δ = 0.025
    γ = 1.0
    φ = 1.6
    α = 0.35
    ρ = 0.9
    A = 10.0 # normalization for ss a
end

ss = get_steady_state(rbc)[:,1]
