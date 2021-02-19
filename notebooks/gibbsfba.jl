using CobraTools
using JuMP
using Gurobi # use your favourite solver
using Measurements
using LinearAlgebra
using JLD
using JSON
using Plots
pyplot()

# E. coli model
modelpath = joinpath("..", "models", "iML1515.json")

model = CobraTools.read_model(modelpath)
gibbs = CobraTools.map_gibbs_rxns(model.rxns) # very slow - rather just import this - will need to reload for other models

ecoli_kJmolCarbon = -37.36 ± 6.20 # formation of biomass kJ/mol upper bound: 74.36 ± 8.67

cbmodel, v, mb, ubs, lbs = CobraTools.CBM(model)
set_optimizer(cbmodel, Gurobi.Optimizer)
set_optimizer_attribute(cbmodel, "OutputFlag", 0) # quiet

biomass_index = model[findfirst(model.rxns, "BIOMASS_Ec_iML1515_core_75p37M")] 
glucose_index = model[findfirst(model.rxns, "EX_glc__D_e")]
o2_index = model[findfirst(model.rxns, "EX_o2_e")]
atpm_index = model[findfirst(model.rxns, "ATPM")]

# Fix glucose use 1.0 then normalization is easy. NB - if not 1 then change normalization!!
CobraTools.set_bound(glucose_index, ubs, lbs; ub=-1.0, lb=-1.0)

# Aerobic
CobraTools.set_bound(o2_index, ubs, lbs; ub=0.0, lb=-1000.0)
# Anaerobic
# CobraTools.set_bound(o2_index, ubs, lbs; ub=0.0, lb=0.0)

# No free ATP generation
CobraTools.set_bound(atpm_index, ubs, lbs; ub=1000.0, lb=0.0)

exts = [(i, model.rxns[i].id)  for i in eachindex(model.rxns) if startswith(model.rxns[i].id, "EX_") && haskey(gibbs, model.rxns[i].id)]

@objective(cbmodel, Min, sum(v[i]*Measurements.value(gibbs[id]) for (i, id) in exts))

optimize!(cbmodel) 
termination_status(cbmodel) != MOI.OPTIMAL && @warn "Optimization issue..."

ΔG_max = objective_value(cbmodel)

fluxes = CobraTools.map_fluxes(v, model)
CobraTools.get_exchanges(fluxes; topN=8, ignorebound=10000)

open(joinpath("escher-fluxes","gibbsfluxes.json"), "w") do io
    JSON.print(io, fluxes)
end