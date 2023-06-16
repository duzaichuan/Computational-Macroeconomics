# using MacroModelling, StatsPlots

@kwdef struct Param
    α = 0.35
    β = 0.9901
    δ = 0.025
    γ = 1.0
    ψ = 1.7333
    ρₐ = 0.9
    η₁ = 2.0
    η₂ = 1.5
end

function rbc_logutil_ss(pa::Param)
    # read-out parameters
    (; α, β, δ, γ, ψ, ρₐ) = pa
    ss = Dict{Symbol, Float64}()

    # compute steady-state analytically
    a = 1.0
    rk = 1.0/β + δ - 1.0
    k_n = ((α * a)/rk)^(1.0/(1.0 - α))
    @assert k_n > 0 "steady-state could not be computed" # impose non-negativity constraints to rule out certain steady-states
    w = (1 - α) * a * k_n^α
    iv_n = δ*k_n
    y_n = a*k_n^α
    c_n = y_n - iv_n
    @assert c_n > 0 "steady-state could not be computed"
    n = γ/ψ * c_n^-1 * w / (1+ γ/ψ * c_n^-1 * w) # closed-form expression for n

    c = c_n * n
    iv = iv_n * n
    k = k_n * n
    y = y_n * n
    uc = γ * c^-1
    un = -ψ/(1-n)
    fn = (1-α) * a * (k/n)^α
    fk = α * a * (k/n)^(α - 1)

    # write to ss dictionary
    ss[:y] = y
    ss[:c] = c
    ss[:k] = k
    ss[:n] = n
    ss[:a] = a
    ss[:rk] = rk
    ss[:w] = w
    ss[:iv] = iv
    ss[:uc] = uc
    ss[:un] = un
    ss[:fn] = fn
    ss[:fk] = fk

    return ss
end

pa = Param()
ss_log = rbc_logutil_ss(pa)

@model rbc_logutil begin
    uc[0] = γ * c[0]^-1
    un[0] = -ψ / (1-n[0])
    fn[0] = (1-α) * α[0] * (k[-1] / n[0])^α
    fk[0] = α * a[0] * (k[-1] / n[0])^(α - 1)

    uc[0] = β * uc[1] * (1 - δ + rk[1])
    w[0] = -un[0]/uc[0]
    w[0] = fn[0]
    rk[0] = fk[0]
    y[0] = a[0] * k[-1]^α * n[0]^(1-α)
    k[0] = (1-δ) * k[-1] + iv[0]
    y[0] = c[0] + iv[0]
    log(a[0]) = ρₐ * log(a[-1]) + eps_a[x]

    # steady_state_model
    a[ss] = 1.0
    rk[ss] = 1/β + δ - 1.0
end

@parameters rbc_logutil begin
    # Targets
    K_Y = 10.0
    IV_Y = 0.25
    WN_Y = 0.65
    N = 1/3

    α = 1 - WN_Y
    δ = IV_Y / K_Y
    RK = α/K_Y
    β = 1/(RK - δ + 1)
    γ = 1.0
    A = 1.0
    KN = ((α * A)/RK)^(1/(1-α))
    K = KN * N
    Y = K/K_Y
    IV = δ * K
    W = (1-α)*A*KN^α
    C = Y - IV
    ψ = γ * (C/N)^-1 * W * (N/(1-N))^-1
    ρₐ = 0.9
end

ss = get_steady_state(rbc_logutil)[:,1]
plot_irf(rbc_logutil, periods=30, variables = [:y, :c, :k, :n, :rk, :w, :iv, :a])

# sims = get_simulation(rbc_logutil, periods=100, variables = [:y, :c, :k, :n, :rk, :w, :iv, :a], initial_state=[1.152, 0.864, 11.517, 0.334, 0.035, 2.246, 0.288, 1.0])

function rbc_ces_ss(pa::Param)
    # read-out parameters
    (; α, β, δ, γ, ψ, ρₐ, η₁, η₂) = pa
    ss = Dict{Symbol, Float64}()

    # compute steady-state analytically
    a = 1.0
    rk = 1.0/β + δ - 1.0
    k_n = ((α * a)/rk)^(1.0/(1.0 - α))
    @assert k_n > 0 "steady-state could not be computed" # impose non-negativity constraints to rule out certain steady-states
    w = (1 - α) * a * k_n^α
    iv_n = δ*k_n
    y_n = a*k_n^α
    c_n = y_n - iv_n
    @assert c_n > 0 "steady-state could not be computed"
    if η₁ == 1.0 && η₂ == 1.0
        n = γ/ψ * c_n^-1 * w / (1+ γ/ψ * c_n^-1 * w) # closed-form expression for n
    else
        f(n) = w * γ* c_n^-η₁ - ψ*(1-n)^-η₂ * n^η₁
        n₀ = 1/3
        n = find_zero(f, n₀)
    end

    # write to ss dictionary
    ss[:y] = y_n * n
    ss[:c] = c_n * n
    ss[:k] = k_n * n
    ss[:n] = n
    ss[:a] = a
    ss[:rk] = rk
    ss[:w] = w
    ss[:iv] = iv_n * n
    ss[:uc] = γ * ss[:c]^-1
    ss[:un] = -ψ/(1-n)
    ss[:fn] = (1-α) * a * (ss[:k]/n)^α
    ss[:fk] = α * a * (ss[:k]/n)^(α - 1)

    return ss
end

ss_ces = rbc_ces_ss(pa)

@model rbc_ces begin
    uc[0] = γ*c[0]^-η₁
    un[0] = -ψ*(1-n[0])^-η₂
    fn[0] = (1-α)*a[0]*(k[-1]/n[0])^α
    fk[0] = α*a[0]*(k[-1]/n[0])^(α-1)
    uc[0] = β*uc[1] *(1-δ+rk[1])
    w[0] = -un[0]/uc[0]
    w[0] = fn[0]
    rk[0] = fk[0]
    y[0] = a[0]*k[-1]^α*n[0]^(1-α)
    k[0] = (1-δ)*k[-1] + iv[0]
    y[0] = c[0] + iv[0]
    log(a[0]) = ρₐ*log(a[-1]) + ϵₐ[x]
end

@parameters rbc_ces begin
    α = 0.35
    β = 0.9901
    δ = 0.025
    γ = 1.0
    ψ = 1.733
    ρₐ = 0.9
    η₁ = 2.0
    η₂ = 1.5
end
get_steady_state(rbc_ces)

@model BK_steady begin
    k[0] = (1-δ_k) * k[-1] + iv[0]
    kg[0] = (1-δ_k) * kg[-1] + ivg[0]
    (1-τ[0]) * w[0] = θ_l * c[0] / (1-n[0])
    c[1]/c[0] = β * (1-δ_k + (1-τ[1])*rk[1])
    y[0] = n[0]^θ_n * k[-1]^θ_k * kg[-1]^θ_g
    w[0]*n[0] = θ_n * y[0]
    rk[0]*k[-1] = θ_k * y[0]
    gb[0] + ivg[0] = τ[0] * (w[0]*n[0] + rk[0]*k[-1]) #- tr[0]
    gb[0] - GB_BAR = ρ_gb*(gb[-1] - GB_BAR) + e_gb[x]
    ivg[0] - IVG_BAR = ρ_ivg * (ivg[-1] - IVG_BAR) + e_ivg[x]
    log(τ[0]/τBAR) = ρ_τ * log(τ[-1]/τBAR) + e_τ[x]
    y[0] = c[0] + iv[0] + gb[0] +ivg[0]
end

@parameters BK_steady begin
    Y_BAR = 1.0 # normalize higher than 1 otherwise consumption might be negative
    GB_BAR = 0.2 * Y_BAR
    IVG_BAR = 0.02 * Y_BAR
    TR_BAR = 0.0
    W_BAR = 2.0
    N_BAR = 1/3
    θ_n = W_BAR*N_BAR/Y_BAR
    θ_k = 1 - θ_n
    θ_g = 0.3 * θ_k
    δ_k = 0.025
    KG_BAR = IVG_BAR/δ_k
    K_BAR = (Y_BAR/(KG_BAR^θ_g * N_BAR^θ_n))^(1/θ_k)
    IV_BAR = K_BAR * δ_k
    C_BAR = Y_BAR - IV_BAR - GB_BAR - IVG_BAR
    RK_BAR = θ_k * Y_BAR/K_BAR
    τBAR = (GB_BAR + IVG_BAR + TR_BAR)/(W_BAR * N_BAR + RK_BAR*K_BAR)
    β = 1/ (1 - δ_k + (1-τBAR)*RK_BAR)
    θ_l = (1- τBAR) * W_BAR * (1-N_BAR)/C_BAR
    ρ_gb = 0.75
    ρ_ivg = 0.75
    ρ_τ = 0.75
    # c > 0
    # w > 0
    # iv > 0
    # rk < 0.05
end

ss_bk = get_steady_state(BK_steady)[:,1]

@model BK_growth begin
    γₓ * k[0] = (1-δₖ) * k[-1] + iv[0]
    γₓ * kg[0] = (1-δₖ) * kg[-1] + ivg[0]
    (1 - τ[0]) * w[0] = θl * c[0] / (1 - n[0])
    γₓ * c[1] / c[0] = β * (1 - δₖ + (1 - τ[1]) * rₖ[1])
    y[0] = n[0]^θₙ * k[-1]^θₖ * kg[-1]^θ_g
    w[0] * n[0] = θₙ * y[0]
    rₖ[0] * k[-1] = θₖ * y[0]
    gb[0] + ivg[0] = τ[0] * (w[0] * n[0] + rₖ[0] * k[-1]) #- tr[0]
    gb[0] - GB_BAR = ρ_gb * (gb[-1] - GB_BAR) + ϵ_gb[x]
    ivg[0] - IVG_BAR = ρ_ivg * (ivg[-1] - IVG_BAR) + ϵ_ivg[x]
    log(τ[0]/τ̄) = ρ_τ * log(τ[-1]/τ̄) + ϵ_τ[x]
    y[0] = c[0] + iv[0] + gb[0] + ivg[0]
end

@parameters BK_growth begin
    γₓ = 1+0.06/4
    Ȳ = 10.0
    GB_BAR = 0.2 * Ȳ
    IVG_BAR = 0.02 * Ȳ
    TR_BAR = 0.0
    W̄ = 2.0
    N̄ = 1/3
    δₖ = 0.025

    KG_BAR = IVG_BAR/(γₓ - 1 + δₖ)
    θₙ = W̄*N̄/Ȳ
    θ_g = 0.3 * (1 - θₙ)
    θₖ = 1 - θₙ - θ_g
    K̄ = (Ȳ/(KG_BAR^θ_g * N̄^θₙ))^(1/θₖ)
    IV_BAR = K̄ * (γₓ -1 + δₖ)
    C̄ = Ȳ - IV_BAR - GB_BAR - IVG_BAR
    RK_BAR = θₖ * Ȳ/K̄
    τ̄ = (GB_BAR + IVG_BAR + TR_BAR)/(W̄ * N̄ + RK_BAR*K̄)

    β = γₓ / (1 - δₖ + (1-τ̄)*RK_BAR)
    θₗ = (1- τ̄) * W̄ * (1-N̄)/C̄
    ρ_gb = 0.75
    ρ_ivg = 0.75
    ρ_τ = 0.75    
end

ss_bk_gr = get_steady_state(BK_growth)[:,1]

@model BK1993_fig2 begin
    uc[0] = c[0]^-1
    ul[0] = θₗ * l[0]^-1
    y[0] = A * k[-1]^θₖ * n[0]^θₙ
    fk[0] = θₖ * A * k[-1]^(θₖ - 1) * n[0]^θₙ
    fn[0] = θₖ * A * k[-1]^θₖ * n[0]^(θₙ - 1)
    γₓ * k[0] = (1- δₖ) * k[-1] + iv[0]
    l[0] + n[0] = 1
    c[0] + iv[0] = (1 - τ[0]) * y[0] + tr[0] + check_walras[0]
    c[0] + iv[0] + gb[0] = y[0]
    τ[0] * y[0] = gb[0] + tr[0]
    uc[0] = λ[0]
    ul[0] = λ[0] * (1 - τ[0]) * fn[0]
    β * λ[1] * (q[1] + 1 - δₖ) = γₓ * λ[0]
    q[0] = (1 - τ[0]) * fk[0]
    gb[0] = GB_BAR + γ_b * eps_gb[x]
    τ[0] = τ̄ + γ_τ * ϵ_τ[x]
    1 + r[0] = γₓ * λ[0] / (λ[1] * β)
    w[0] = fn[0]
end


@parameters BK1993_fig2 begin
    A = 1.0
    γ_x = 1.016
    γ_b = 1.0
    γ_τ = 0.0
    θ_k = 0.42
    θ_n = 1 - θ_k
    δ_k = 0.1
    N = 0.2
    L = 1 - N
    R = 1.065
    β = γ_x * R^-1
    sG = 0.2
    τBAR = 0.2
    τ = τBAR
    Q = γ_x / β -1 + δ_k
    FK = Q/(1 - τBAR)
    K = (FK / ((θ_k * A * N^θ_n)))^(1 / (θ_k - 1))
    W = θ_n * A * K^θ_k * N^(θ_n - 1)
    IV = (γ_x - 1 + δ_k) * K
    Y = A * N^(1- θ_k) * K^θ_k
    GB_BAR = sG * Y
    C = Y - IV - GB_BAR
    UC = C^-1
    UL = UC * (1 - τBAR) * W
    θ_l = UL * L
    τ[ss] -> τBAR
end

ss_bk_fig2 = get_steady_state(BK1993_fig2)[:,1]

