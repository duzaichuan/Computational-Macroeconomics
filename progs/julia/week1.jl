using MacroModelling, StatsPlots
# Random.seed!(123)
# !! seed set does not work for get_simulation

@model arma11 begin
    x[0] - θ*x[-1] = ϵ[exo] - φ*ϵ[exo-1]
end

@parameters arma11 begin
    φ = 0.4
    θ = 0.4
end

plot_simulations(arma11, periods=200, initial_state=Float64.(arma11.solved_vals))

sims = get_simulation(arma11, periods=200, initial_state=Float64.(arma11.solved_vals))

dyn_x = sims[1,:]
dyn_e = sims[2,:]
sample_size = length(sims.Periods)

θ = arma11.parameter_values[2]
φ = arma11.parameter_values[1]
histogram(dyn_e)

e = randn(sample_size)
if isequal(e, dyn_e)
    println("shock series are equal")
else
    println("shock series are not equal!")
end

function arma_sim(e)
    x = zeros(length(e)) # preallocate storage
    x[1] = θ*0 + e[1] - φ*0
    for t = 2:length(e)
        x[t] = θ * x[t-1] + e[t] - φ * e[t-1]
    end
    x
end

x = arma_sim(e)
if isequal(x, dyn_x)
    println("simulated series are equal")
else
    println("simulated series are not equal!")
end
