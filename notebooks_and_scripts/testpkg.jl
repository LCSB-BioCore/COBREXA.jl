using CobraTools
using JuMP
using Tulip

model = CobraTools.read_model(joinpath("models", "e_coli_core.json"))
# optimizer = Tulip.Optimizer
# biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
# sol = fba(model, biomass, Tulip.Optimizer) # classic flux balance analysis

# atts = Dict("OutputFlag" => 0)
# fva_max, fva_min = fva(model, biomass, Gurobi.Optimizer, solver_attributes=atts)
# fva_max, fva_min = fva(model, biomass, Tulip.Optimizer)
# optimizer = Tulip.Optimizer 
# atts = Dict("IPM_IterationsLimit" => 400)
# fva_max, fva_min = fva(model, biomass, optimizer, solver_attributes=atts)
# sol = fba(model, biomass, Tulip.Optimizer) # classic flux balance analysis
# atom_exchange(sol, model)
# consuming, producing = exchange_reactions(sol; verbose=false)
# consuming, producing = metabolite_fluxes(sol, model)

##################################

# model = CobraTools.read_model(joinpath("models", "e_coli_core.json"))
# optimizer = Tulip.Optimizer
# biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
# cons = Dict("EX_glc__D_e" => (-12.0, -12.0))
# sol = fba(model, biomass, optimizer, constraints=cons) # classic flux balance analysis
# cons["BIOMASS_Ecoli_core_w_GAM"] = (sol["BIOMASS_Ecoli_core_w_GAM"], sol["BIOMASS_Ecoli_core_w_GAM"]*0.99)

# # sample!
# (@view shat[:])

# samples = @time CobraTools.hit_and_run(100_000, model, optimizer; keepevery=10, samplesize=5000, constraints=cons)
# println(mean(samples[64,:]))
# println(std(samples[64,:]))

# samples = @time CobraTools.achr(100_000, model, optimizer; keepevery=10, samplesize=5000, constraints=cons)
# println(mean(samples[64,:]))
# println(std(samples[64,:]))


brenda_data = parse_brenda(joinpath("test", "data", "small_brenda.txt"))


# brenda_data = parse_brenda(joinpath("data", "small_brenda.txt"))
# model = CobraTools.read_model(joinpath("models", "e_coli_core.json"))
# biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
# atts = Dict("eps_abs" => 5e-4,"eps_rel" => 5e-4, "max_iter" => 100_000, "verbose"=>false) 
# sol = pfba(model, biomass, [Tulip.Optimizer, OSQP.Optimizer]; solver_attributes=Dict("opt1" => Dict{Any, Any}(), "opt2" => atts)) # try two optimizers
 
# # sol = pfba(model, biomass, [Tulip.Optimizer, ]) # classic flux balance analysis

# rs = CobraTools.build_rxn_string(model.reactions[64]) == "1.0 KEGG:C00354 = 1.0 KEGG:C00111 + 1.0 KEGG:C00661"


# gibbs_str = JSON.parsefile(joinpath("test", "data", "gibbs.json"))
# gibbs = Dict{String, Measurement{Float64}}()
# for (k, v) in gibbs_str
#     gibbs[k] = parse(Measurement{Float64}, v)
# end

# total_ΔG, mf = map_gibbs_external(sol, gibbs)

# total_ΔG, mf = map_gibbs_internal(sol, gibbs) 

# for (k, v) in gibbs
#     if abs(sol[k]) > 0.000 
#         println(k, "\t", round(v*sol[k], digits=3))
#     end
# end
