using MacroModelling, StatsPlots
# import SymPy as sym
using Symbolics, SymbolicNumericIntegration

# 1
1//3

# 2
sin(π)

# 3
@variables x y z

# 4

# 6
@variables a b c # Type Num, subtype of Real (eg. 1//3)
# @syms a b c
A = [a b c; c a b; b c a]
isequal(sum(A[1,:]), sum(A[:, 2]))

# 8
@variables a b c x
g_l(x) = a*x + b*x +c
Symbolics.solve_for(g_l(x), x) # Symbolics only feasible for linear solve
# @syms a b c x
# solve(g(x),x)
# solve(f-2,x)
f = a*x^2 + b*x + c
g(x) = a*x^2 + b*x +c
g(4)

# 9
Symbolics.gradient(f, [a,b,c,x])
Symbolics.derivative(f, a)
Symbolics.derivative(f, x)

# 10
int = integrate(f, x)[1]
res = substitute(int, Dict(x => 1)) - substitute(int, Dict(x => -2))
simplify(res, expand = true)
# or
f_int_ex = build_function(int, x)
f_int = eval(f_int_ex)
f_int(1) - f_int(-2)

# 11
@variables x
φ = (1 + sqrt(x))/2
golden_ratio = φ^2 - φ - 1
simplify(golden_ratio, expand = true)

# 12
g = (x^2−1)*(x^4+x^3+x^2+x+1)*(x^4−x^3+x^2−x+1)
simplify(g, expand = true)

# 13

# 14


@model NK begin
    rnom[0] = rreal[0] * pie[1]

    # households
    k[0] = (1 - δ) * k[-1] + (1 - φ_iv/2 * (iv[0]/iv[-1] - 1)^2) * iv[0]
    λ[0] = ζ[0] * c[0]^(-σ_c)
    w[0] = χ_h * h[0]^σ_h * c[0]^σ_c
    λ[0] = β * λ[1] * rreal[0]
    1 = q[0] - q[0] * φ_iv * (iv[0]/iv[-1] - 1) * (iv[0]/iv[-1]) - q[0] * φ_iv/2 * (iv[0]/iv[-1] - 1)^2 + β * λ[1]/ λ[0] * q[1] * φ_iv * (iv[1]/iv[0] - 1) * (iv[1]/iv[0])^2
    q[0] = β * λ[1] / λ[0] * (rk[1] + q[1] * (1 - δ))

    # firms
    k[-1] / h[0] = w[0] / (1 - α) * (α/rk[0])
    mc[0] = 1 / a[0] * (w[0]/ (1 - α))^(1-α) * (rk[0]/α)^α
    ptilde[0] * s1[0] = ϵ / (ϵ - 1) * s2[0]
    s1[0] = y[0] + β * θ * λ[1] / λ[0] * pie[1]^(ϵ - 1) * s1[1]
    s2[0] = mc[0] * y[0] + β * θ * λ[1] / λ[0] * pie[1]^ϵ * s2[1]
    div[0]= y[0] - w[0]*h[0] - rk[0] * k[-1]
    1 = θ * pie[0]^(ϵ - 1) + (1-θ) * ptilde[0]^(1-ϵ)

    # aggregation and market clear
    y[0] = c[0] + iv[0]
    pstar[0] * y[0] = a[0] * k[-1]^α * h[0]^(1-α)
    pstar[0] = (1-θ) * ptilde[0]^(-ϵ) + θ * pie[0]^ϵ * pstar[-1]
    rnom[0] = rnom[ss] * (pie[0]/PISTAR)^ψ_pi * (y[0]/y[ss])^ψ_y * exp(ν[0])

    # exogenous processes
    log(a[0]) = ρ_a * log(a[-1]) + eps_a[x]
    ν[0] = ρ_ν * ν[-1] + eps_ν[x]
    log(ζ[0]) = ρ_ζ * log(ζ[-1]) + eps_ζ[x]

    # reporting
    hat_y[0] = log(y[0]) - log(y[ss])
    hat_w[0] = log(w[0]) - log(w[ss])
    hat_h[0] = log(h[0]) - log(h[ss])
    hat_pi_ann[0] = 4*(log(pie[0]) - log(pie[ss]))
    hat_rnom_ann[0] = 4*(log(rnom[0]) - log(rnom[ss]))
    hat_rreal_ann[0] = 4*(log(rreal[0]) - log(rreal[ss]))
    hat_mc[0] = log(mc[0]) - log(mc[ss])
    hat_a[0] = log(a[0]) - log(a[ss])
    hat_ζ[0] = log(ζ[0]) - log(ζ[ss])
end

@parameters NK begin
    δ = 0.025
    φ_iv = 4.25
    σ_c = 2.0
    σ_h = 2.0
    β = 0.99
    α = 0.33
    θ = 0.75
    ϵ = 6.0
    χ_h = 5.0
    PISTAR = 1.005
    ψ_pi = 1.5
    ψ_y = 0.5/4
    ρ_ζ = 0.5
    ρ_a = 0.9
    ρ_ν = 0.5
end

ss_nk = get_steady_state(NK)[:,1]
