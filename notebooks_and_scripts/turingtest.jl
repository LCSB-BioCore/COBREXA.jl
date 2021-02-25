using CobraTools
using JuMP
using Gurobi # use your favourite solver
using Measurements
using LinearAlgebra
using MCMCChains
using HypothesisTests
using StatsBase
using Turing
using Distributions
using DataFrames
using StatsPlots
pyplot()

# E. coli model
modelpath = joinpath("models", "e_coli_core.json") 
model = CobraTools.read_model(modelpath)

#### set bounds
cbmodel, v, mb, ubs, lbs = CobraTools.CBM(model)
set_optimizer(cbmodel, Gurobi.Optimizer)
set_optimizer_attribute(cbmodel, "OutputFlag", 0) # quiet

biomass_index = model[findfirst(model.rxns, "BIOMASS_Ecoli_core_w_GAM")] 
glucose_index = model[findfirst(model.rxns, "EX_glc__D_e")]
o2_index = model[findfirst(model.rxns, "EX_o2_e")]
atpm_index = model[findfirst(model.rxns, "ATPM")]
etoh_index = model[findfirst(model.rxns, "EX_etoh_e")]

# Fix glucose use 1.0 then normalization is easy. NB - if not 1 then change normalization!!
CobraTools.set_bound(glucose_index, ubs, lbs; ub=-1.0, lb=-1.01)

# Aerobic
# CobraTools.set_bound(o2_index, ubs, lbs; ub=1000.0, lb=-1000.0)
# Anaerobic
CobraTools.set_bound(o2_index, ubs, lbs; ub=1000.0, lb=0.0)

# No free ATP generation
CobraTools.set_bound(atpm_index, ubs, lbs; ub=1000.0, lb=0.0)

@objective(cbmodel, Max, v[biomass_index])

optimize!(cbmodel) 
termination_status(cbmodel) != MOI.OPTIMAL && @warn "Optimization issue..."

μ = objective_value(cbmodel)

### Fix biomass as a constraint
CobraTools.set_bound(biomass_index, ubs, lbs; ub=μ, lb=0.99*μ)


ubvec, lbvec = CobraTools.get_bound_vectors(ubs, lbs)
S, _, _, _ = CobraTools.get_core_model(model)

###########################################
#### Do inference

S = [1 -1.0 0 0 0 0 0;0 1 -1 -1 -1 0 0;0 0 1 1 1 -1 0;0 0 0 0 0 1 -1]
Sin = S[:, 2:end]
@model function turingtest(S)

    Σ ~ InverseWishart(6 + 1, Matrix{Float64}(I, 6, 6))    
    
    flux ~ MvNormal(zeros(6), Symmetric(Σ))

    θ ~ Uniform(0.0, 100.0)
    b = sum(abs, S*flux - [1.0, 0.0, 0.0, 0.0]) # input is r₀ = -1
    0.0 ~ truncated(Normal(b, θ), 0.0, Inf)
end

tmodel = turingtest(Sin)
chain = sample(tmodel, NUTS(0.65), 1_000)

################

df = DataFrame(x1 = rand(Normal(), 100), x2 = rand(Normal(), 100), x3 = rand(Normal(), 100))
beta = [5, 4, 5, 6]
df[!, :y] = beta[1] .+ beta[2]*df[!, :x1] + beta[3]*df[!, :x2] + beta[4]*df[!, :x3] + rand(Normal(), 100)

@model BayesReg(y, x1, x2, x3) = begin
    
    n_obs = length(y)

    σ₂ ~ InverseGamma(0.001, 0.001)

    Sigma   ~ InverseWishart(
        4 + 1, 
        Matrix{Float64}(I, 4, 4) 
    )
    Mu ~ MvNormal( 
        zeros(4),    
        Matrix{Float64}(I, 4, 4)   
    )

    beta_hat ~ MvNormal(Mu, Symmetric(Sigma))

    mu = beta_hat[1] .+ beta_hat[2]*x1 .+ beta_hat[3]*x2 .+ beta_hat[4]*x3 
    for i = 1:n_obs
        y[i] ~ Normal(mu[i], σ₂)
    end
end


model = BayesReg(df[!, :y], df[!, :x1], df[!, :x2], df[!, :x3])
chain = sample(model, NUTS(0.65), 1000)

