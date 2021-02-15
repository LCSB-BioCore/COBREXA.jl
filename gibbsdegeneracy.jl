using CobraTools
using JuMP
using Gurobi # use your favourite solver
using Measurements
using LinearAlgebra
using JLD
using Plots
pyplot()

# E. coli model
modelpath = joinpath("models", "iJO1366.json") 
model = CobraTools.read_model(modelpath)

decomp = JLD.load(joinpath("data", "dgzeros.jld"), "gibbs")
gibbs = Dict{String, Measurement{Float64}}()
for (k, vs) in decomp
    gibbs[k] = vs[1] ± vs[2]
end
##

ecoli_kJmolCarbon = -37.36 ± 8.55 # formation of biomass kJ/mol 74.36 ± 8.67

cbmodel, v, mb, ubs, lbs = CobraTools.CBM(model)
set_optimizer(cbmodel, Gurobi.Optimizer)
set_optimizer_attribute(cbmodel, "OutputFlag", 0) # quiet

biomass_index = model[findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")] 
glucose_index = model[findfirst(model.rxns, "EX_glc__D_e")]
o2_index = model[findfirst(model.rxns, "EX_o2_e")]
atpm_index = model[findfirst(model.rxns, "ATPM")]

# Fix glucose use 1.0 then normalization is easy. NB - if not 1 then change normalization!!
CobraTools.set_bound(glucose_index, ubs, lbs; ub=-1.0, lb=-1.0)

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

wpoints = CobraTools.get_warmup_points(cbmodel, v, ubs, lbs) # very slow

Nsamples = 100_000
samples = CobraTools.achr(Nsamples, wpoints, model, ubs, lbs, burnin=0.1, keepevery=200) 

ΔG_exts = Measurement{Float64}[]
for s in 1:size(samples, 2)
    fluxes = CobraTools.map_fluxes(samples[:, s], model)
    carbon_ex = CobraTools.atom_exchange(fluxes, model)["C"] # carbon flux
    ΔG_ext, missing_ext =  CobraTools.map_gibbs_external(fluxes, gibbs) 
    ΔG_ext -= carbon_ex*ecoli_kJmolCarbon # minus because carbons consumed
    push!(ΔG_exts, ΔG_ext)
end