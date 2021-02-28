using CobraTools
using JuMP
using Gurobi
using Tulip

model = CobraTools.read_model(joinpath("models", "e_coli_core.json"))
biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
# sol = fba(model, biomass, Tulip.Optimizer) # classic flux balance analysis
# atom_exchange(sol, model)

# consuming, producing = exchange_reactions(sol; verbose=false)

# consuming, producing = metabolite_fluxes(sol, model)

##################################
# Get warmup points


cbmodel, v, mb, ubs, lbs = CobraTools.build_cbm(model)
set_optimizer(cbmodel, Tulip.Optimizer)

biomass_index = model[findfirst(model.rxns, "BIOMASS_Ecoli_core_w_GAM")] 
glucose_index = model[findfirst(model.rxns, "EX_glc__D_e")]
o2_index = model[findfirst(model.rxns, "EX_o2_e")]
atpm_index = model[findfirst(model.rxns, "ATPM")]

CobraTools.set_bound(glucose_index, ubs, lbs; ub=-1.0, lb=-1.0)
CobraTools.set_bound(o2_index, ubs, lbs; ub=1000.0, lb=0.0)
CobraTools.set_bound(atpm_index, ubs, lbs; ub=1000.0, lb=0.0)

# Do FBA
@objective(cbmodel, Max, v[biomass_index])
optimize!(cbmodel) 
μ_max = round(objective_value(cbmodel), digits=6)

CobraTools.set_bound(biomass_index, ubs, lbs; ub=μ_max, lb=0.9*μ_max)


wpoints = CobraTools.get_warmup_points(cbmodel, v, ubs, lbs, numstop=4) # very slow

# sample!
# samples = @time CobraTools.hit_and_run(100_000, wpoints, ubs, lbs; keepevery=10, samplesize=5000, W=1000) 

###########################
# violation_inds = CobraTools.test_samples(samples, model, ubs, lbs)

