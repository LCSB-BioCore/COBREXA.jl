using CobraTools
using Gurobi

using JuMP
using Tulip


model = Model(Tulip.Optimizer)
@variable(model, -1.0 <= x[1:2] <= 2)
@objective(model, Max, 5x[1] + 3 * x[2])
@constraint(model, con1, 1x[1] + 5x[2] <= 3)
@constraint(model, con2, 0.0 .<= x)


optimize!(model)
objective_value(model)

set_normalized_rhs(con2[1], 0.2)
set_normalized_rhs(con2[2], 0.2)


model = CobraTools.read_model(joinpath("models", "e_coli_core.json"))
biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
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
model = CobraTools.read_model(joinpath("data", "iJO1366.json") )
    
    rxn_original = findfirst(model.reactions, "NADH16pp")
    nadh = findfirst(model.metabolites, "nadh_c")
    h_c = findfirst(model.metabolites, "h_c")
    q8 = findfirst(model.metabolites, "nadp_c")
    q8h2 = findfirst(model.metabolites, "g3p_c")
    nad = findfirst(model.metabolites, "nad_c")
    h_p = findfirst(model.metabolites, "for_c")
    
    rxn = nadh + 4.0*h_c + 1.0*q8 ⟶  1.0*q8h2 + 1.0*nad + 3.0*h_p
    @test rxn.lb == 0.0 && rxn.ub > 0.0

    rxn = 1.0*nadh + 4.0*h_c + q8 ← 1.0*q8h2 + 1.0*nad + 3.0*h_p
    @test rxn.lb < 0.0 && rxn.ub == 0.0
  
    rxn = 1.0*nadh + 4.0*h_c + 1.0*q8 ↔ q8h2 + nad + 3.0*h_p
    @test rxn.lb < 0.0 && rxn.ub > 0.0
    
    rxn = 1.0*nadh → ∅
    @test length(rxn.metabolites) == 1

    rxn = ∅ → nadh
    @test length(rxn.metabolites) == 1

    rxn = ∅ → 1.0nadh
    @test length(rxn.metabolites) == 1
    
    rxn = 1.0*nadh + 4.0*h_c + 1.0*q8 ⟶  1.0*q8h2 + 1.0*nad + 3.0*h_p
    @test prod(values(rxn.metabolites)) == -12
    @test ("q8h2_c" in [x.id for x  in keys(rxn.metabolites)])

    rxn = nadh + 4.0*h_c + 1.0*q8 ⟶  1.0*q8h2 + 1.0*nad + 3.0*h_p
    @test rxn.lb == 0.0 && rxn.ub > 0.0