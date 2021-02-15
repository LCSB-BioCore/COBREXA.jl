using CobraTools
using Measurements
using LinearAlgebra
using Statistics
using JLD
using JuMP
using Gurobi

# E. coli model
modelpath = joinpath("models", "iJO1366.json") 
model = CobraTools.read_model(modelpath)

###### Set constraints
cbmodel, v, mb, ubs, lbs = CobraTools.CBM(model)
set_optimizer(cbmodel, Gurobi.Optimizer)
set_optimizer_attribute(cbmodel, "OutputFlag", 0) # quiet

biomass_index = model[findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")] 
glucose_index = model[findfirst(model.rxns, "EX_glc__D_e")]
o2_index = model[findfirst(model.rxns, "EX_o2_e")]
atpm_index = model[findfirst(model.rxns, "ATPM")]
etoh_index = model[findfirst(model.rxns, "EX_etoh_e")]


CobraTools.set_bound(glucose_index, ubs, lbs; ub=-1.0, lb=-1.0)
CobraTools.set_bound(o2_index, ubs, lbs; ub=1000.0, lb=0.0)
CobraTools.set_bound(atpm_index, ubs, lbs; ub=1000.0, lb=0.0)

# Do FBA
@objective(cbmodel, Max, v[biomass_index])
optimize!(cbmodel) 
termination_status(cbmodel) != MOI.OPTIMAL && @warn "Optimization issue..."
μ_max = round(objective_value(cbmodel), digits=6)

CobraTools.set_bound(biomass_index, ubs, lbs; ub=μ_max, lb=0.9*μ_max)

##################################
# Get warmup points
wpoints = CobraTools.get_warmup_points(cbmodel, v, ubs, lbs) # very slow

# sample!
samples = @time CobraTools.achr(100_000, wpoints, model, ubs, lbs) 

###########################










