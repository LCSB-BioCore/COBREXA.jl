using CobraTools
using Gurobi

using SBML

m = readSBML(joinpath("models", "iJO1366.xml"))



using CobraTools # hide
met1 = Metabolite()
met1.id = "met1"
met1.name = "Metabolite 1"
met1.formula = "C6H12O6N"
met1.charge = 1
met1.compartment = "c"
met1.notes = Dict("notes"=>["This is a made up metabolite", "Another note"])
met1.annotation = Dict("sboterm" => "sbo000001", "kegg.compound" => ["C0001", "C0010"])

met2 = Metabolite("met2")
met2.formula = "C6H12O6N"

met3 = Metabolite("met3")
met3.formula = "X"
met3.annotation = Dict("sboterm" => "sbo00001", "kegg.compound" => ["C02222", "C0001"])

mets = [met1, met2, met3]

dup, ind = check_duplicate_annotations(mets, met3)
if dup
    println("Duplicate found at index: ", ind)
end

mms = check_same_formula([met3, met1], met2)
println("Metabolites with the same formula as \"met2\":")
mms[1]







# using JuMP
# using Tulip

# 1+1

# model = CobraTools.read_model(joinpath("models", "e_coli_core.json"))
# biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
# sol = fba(model, biomass, Tulip.Optimizer) # classic flux balance analysis
# atom_exchange(sol, model)

# consuming, producing = exchange_reactions(sol; verbose=false)

# consuming, producing = metabolite_fluxes(sol, model)

##################################
# Get warmup points


# cbmodel, v, mb, ubs, lbs = CobraTools.build_cbm(model)
# set_optimizer(cbmodel, Tulip.Optimizer)

# biomass_index = model[findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")] 
# glucose_index = model[findfirst(model.reactions, "EX_glc__D_e")]
# o2_index = model[findfirst(model.reactions, "EX_o2_e")]
# atpm_index = model[findfirst(model.reactions, "ATPM")]

# set_bound(glucose_index, ubs, lbs; ub=-1.0, lb=-1.0)
# set_bound(o2_index, ubs, lbs; ub=1000.0, lb=0.0)
# set_bound(atpm_index, ubs, lbs; ub=1000.0, lb=0.0)

# Do FBA
# @objective(cbmodel, Max, v[biomass_index])
# optimize!(cbmodel) 
# μ_max = round(objective_value(cbmodel), digits=6)

# set_normalized_rhs(ubs[25], 0.1)
# set_bound(biomass_index, ubs, lbs; ub=μ_max, lb=0.9*μ_max)


# wpoints = CobraTools.get_warmup_points(cbmodel, v, ubs, lbs, numstop=4) # very slow

# sample!
# samples = @time CobraTools.hit_and_run(100_000, wpoints, ubs, lbs; keepevery=10, samplesize=5000, W=1000) 

###########################
# violation_inds = CobraTools.test_samples(samples, model, ubs, lbs)


##################3