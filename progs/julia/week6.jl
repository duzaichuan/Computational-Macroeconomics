using MacroModelling, StatsPlots

# predetermined_variables k i s r dd pop;

@model NK_sir begin
    # [name='1) aggregate supply']
    y[0] = pbreve[0] * A * k[-1]^(1-α) * n[0]^α
    # rk[0] - δ = pbreve[0] * A * (1-α) * k[-1]^-α * n[0]^α
    # [name='2) marginal costs']
    mc[0] = w[0]^α * rk[0]^(1-α) / ( A * α^α * (1-α)^(1-α) )
    # [name='3) optimal factor input']
    w[0] = mc[0] * α * A * n[0]^(α-1) * k[-1]^(1-α)
    # [name='4) capital accumulation']
    k[0] = (1-δ) * k[-1] + x[0]
    # [name='5) aggregate demand']
    y[0] = c[0] + x[0] + g_ss
    # [name='6) aggregate hours']
    n[0] = (s[-1] * ns[0]) + (i[-1] * ni[0]) + (r[-1] * nr[0])
    # [name='7) aggregate consumption']
    c[0] = (s[-1] * cs[0]) + (i[-1] * ci[0]) + (r[-1] * cr[0])
    # [name='8) transmission function for new infections']
    τ[0] = π1 * (s[-1] * cs[0]) * (i[-1] * ci[0]) + π2 * (s[-1] * ns[0]) * (i[-1] * ni[0]) + π3 * (s[-1] * i[-1])
    # [name='9) law of motion susceptible population']
    s[0] = s[-1] - τ[0]
    # [name='10) law of motion infected population']
    i[0] = i[-1] + τ[0] - (π_r + π_d) * i[-1]
    # [name='11) law of motion recovered population']
    r[0] = r[-1] + π_r * i[-1]
    # [name='12) law of motion deceased population']
    dd[0] = dd[-1] + π_d * i[-1]
    # [name='13) total population']
    pop[0] = pop[-1] - π_d * i[-1]
    # [name='14) first order condition wrt consumption susceptibles']
    1 / cs[0] = lambtilde[0] - lamtau[0] * π1 * (i[-1] * ci[0])
    # [name='15) first order condition wrt consumption infected']
    1 / ci[0] = lambtilde[0]
    # [name='16) first order condition wrt consumption recovered']
    1 / cr[0] = lambtilde[0]
    # [name='17) first order condition wrt hours susceptibles']
    θ * ns[0] = lambtilde[0] * w[0] + lamtau[0] * π2 * (i[-1] * ni[0])
    # [name='18) first order condition wrt hours infected']
    θ * ni[0] = lambtilde[0] * w[0]
    # [name='19) first order condition wrt hours recovered']
    θ * nr[0] = lambtilde[0] * w[0]
    # [name='20) first order condition wrt capital']
    lambtilde[0] = β * ( rk[+1] + 1 - δ ) * lambtilde[+1]
    # [name='21) first order condition wrt newly infected']
    lami[0] = lamtau[0] + lams[0]
    # [name='22) first order condition wrt susceptibles']
    log(cs[+1]) - θ/2 * ns[+1]^2 + lamtau[+1] * ( π1 * cs[+1] * i[0] * ci[+1] + π2 * ns[+1] * i[0] * ni[+1] + π3 * i[0]) + lambtilde[+1] * ( w[+1] * ns[+1] - cs[+1]) - lams[0] / β + lams[+1]
    # [name='23) first order condition wrt infected']
    log(ci[+1]) - θ/2 * ni[+1]^2 + lambtilde[+1] * ( w[+1] * ni[+1] - ci[+1] ) - lami[0]/β + lami[+1] * (1 - π_r - π_d) + lamr[+1]*π_r = 0
    # [name='24) first order condition wrt recovered']
    log(cr[+1]) - θ/2 * nr[+1]^2 + lambtilde[+1] * ( w[+1] * nr[+1] - cr[+1] ) - lamr[0]/β + lamr[+1] = 0
    # [name='25) first order condition wrt bonds']
    lambtilde[0] = β * rr[0] * lambtilde[+1]
    # [name='26) real interest rate']
    rr[0] = Rb[0] / π_e[+1]
    # [name='27) recursion nonlinear price setting 1']
    Kf[0] = γ * mc[0] * lambtilde[0] * y[0] + β * ξ * π_e[+1]^(γ/(γ-1)) * Kf[+1]
    # [name='28) recursion nonlinear price setting 2']
    F[0] = lambtilde[0] * y[0] + β * ξ * π_e[+1]^(1/(γ-1)) * F[+1]
    # [name='29) nonlinear price setting']
    Kf[0] = F[0] * ( ( 1 - ξ * π_e[0]^(1/(γ-1)) ) / (1 - ξ) )^(-(γ-1))
    # [name='30) inverse price dispersion']
    pbreve[0]^-1 = (1-ξ) * ( ( 1 - ξ * π_e[0]^(1/(γ-1)) ) / (1 - ξ) )^γ + ξ * π_e[0]^(γ/(γ-1)) / pbreve[-1]
    # [name='31) Taylor rule']
    Rb[0] - Rb_ss = r_π * log( π_e[0] / π_e_ss) + r_x * log( y[0] / y_ss)
end


@parameters NK_sir begin
    ξ             = 0.98
    r_π           = 1.5
    r_x           = 0.5/52
    γ             = 1.35
    π1            = 2.568e-7
    π2            = 1.593e-4
    π3            = 0.4997
    π_d           = 7*0.002/14
    π_r           = 7/14 - π_d
    β             = 0.98^(1/52)
    α             = 2/3
    δ             = 0.06/52
    inc_target    = 58000/52
    n_target      = 28.0
    η             = 0.19
    ε             = 0.001
    RplusD_target = 0.6
    π_e_ss        = 1.0
    y_ss          = inc_target
    Rb_ss         = π_e_ss/β
    g_ss          = η*y_ss
    n_ss          = n_target
    mc_ss         = 1/γ
    rk_ss         = 1/β-1+δ
    w_ss          = mc_ss*α*y_ss/n_ss
    kn_ss         = (1-α)*w_ss/α/rk_ss
    yk_ss         = (y_ss/n_ss)/kn_ss
    A             = (y_ss/n_ss)^α*yk_ss^(1-α)
    ns_ss         = n_ss
    k_ss          = (y_ss/A/n_ss^α)^(1/(1-α))
    x_ss          = δ*k_ss
    c_ss          = y_ss - x_ss - g_ss
    cs_ss         = c_ss
    lambtilde_ss  = 1/cs_ss
    θ             = lambtilde_ss*w_ss/ns_ss
    0.0015 < rk < 0.00155
    0.0014 < lambtilde < 0.00144
    0.73 < mc < 0.75
    1.00038 < Rb < 1.00040
    1.00038 < rr < 1.00040
    -0.0001 < π_e < 0.00001
    178550 < k < 178552
    27 < n < 29
    19 < w < 20
    1115 < y < 11116
    206 < x < 207
    697 < c < 698
    697 < cs < 698
    697 < ci < 698
    697 < cr < 698
    0.99 < s < 1.000001
    -0.00001 < i < 0.00001
    -0.00001 < r < 0.00001
    27 < ns < 29
    27 < ni < 29
    27 < nr < 29
    -31 < lamtau < -30
    15261 < lams < 15293
    15261 < lami < 15262
    15261 < lamr < 15293
    0.99 < pbreve < 1.001
    78 < Kf < 79
    78 < F < 79
    # y[ss] -> y_ss
    # Rb[ss] -> Rb_ss
    # π_e[ss] -> π_e_ss
    # n[ss] -> n_ss
    # x[ss] -> x_ss
    # mc[ss] -> mc_ss
    # pop[ss] -> 1.0
    # s[ss] -> 1.0
    # i[ss] -> 0.0
    # r[ss] -> 0.0
    # pbreve[ss] -> 1.0
    # cs[ss] -> cs_ss
    # nr[ss] -> ns_ss
    # lambtilde[ss] -> lambtilde_ss
    # Kf[ss] -> 1/(1−β*ξ)*γ*mc_ss*lambtilde_ss*y_ss
    # F[ss] -> 1/(1−β*ξ)*lambtilde_ss*y_ss
    # τ[ss] -> 0.0
    # n[ss] = n_ss | θ
    # y[ss] = y_ss | A
end
