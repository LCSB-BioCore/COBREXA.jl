using CobraTools
using JuMP
using Gurobi # use your favourite solver
using Measurements
using LinearAlgebra
using MCMCChains
using StatsBase
using Plots
pyplot()

# E. coli model
modelpath = joinpath("models", "iML1515.json") 
model = CobraTools.read_model(modelpath)
gibbs = CobraTools.map_gibbs_rxns(model.rxns) # very slow - rather just import this - will need to reload for other models


ecoli_kJmolCarbon = -37.36 ± 8.55 # formation of biomass kJ/mol 74.36 ± 8.67

cbmodel, v, mb, ubs, lbs = CobraTools.CBM(model)
set_optimizer(cbmodel, Gurobi.Optimizer)
set_optimizer_attribute(cbmodel, "OutputFlag", 0) # quiet

biomass_index = model[findfirst(model.rxns, "BIOMASS_Ec_iML1515_core_75p37M")] 
glucose_index = model[findfirst(model.rxns, "EX_glc__D_e")]
o2_index = model[findfirst(model.rxns, "EX_o2_e")]
atpm_index = model[findfirst(model.rxns, "ATPM")]
etoh_index = model[findfirst(model.rxns, "EX_etoh_e")]

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

##################################
# Get warmup points
wpoints = CobraTools.get_warmup_points(cbmodel, v, ubs, lbs, numstop=1000) # very slow

# sample
samples = @time CobraTools.hit_and_run(1000_000, wpoints, ubs, lbs; keepevery=10, samplesize=5000) 
samples = @time CobraTools.achr(100_000, wpoints, ubs, lbs; keepevery=10, samplesize=5000) 

etoh_chain = Chains(samples[etoh_index, :])
gewekediag(etoh_chain)[1]

plot(samples[etoh_index, :])
v = StatsBase.autocor(samples[etoh_index, :])
plot(v)
###########################
violation_inds = CobraTools.test_samples(samples, model, ubs, lbs)
###########################

ΔG_exts = Measurement{Float64}[]
for s in 1:size(samples, 2)
    fluxes = CobraTools.map_fluxes(samples[:, s], model)
    carbon_ex = CobraTools.atom_exchange(fluxes, model)["C"] # carbon flux
    ΔG_ext, missing_ext =  CobraTools.map_gibbs_external(fluxes, gibbs) 
    ΔG_ext -= carbon_ex*ecoli_kJmolCarbon # minus because carbons consumed
    push!(ΔG_exts, ΔG_ext) # for quick histogram
end

dgs = [x.val for x in ΔG_exts]
histogram(dgs)