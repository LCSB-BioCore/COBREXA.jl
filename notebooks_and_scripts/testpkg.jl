using CobraTools
using JuMP
using Gurobi
using Tulip
using OSQP
using GLPK

# model = read_model(joinpath("models", "e_coli_core.json"))
# biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")

# # optimizer = Gurobi.Optimizer
# # atts = Dict("OutputFlag" => 0)
# # sol = fba(model, biomass, optimizer; solver_attributes=atts)
# # sol["BIOMASS_Ecoli_core_w_GAM"]
# # sol = pfba(model, biomass, optimizer; solver_attributes=atts)


# optimizer = Tulip.Optimizer # quiet by default
# sol = fba(model, biomass, optimizer)
# sol["BIOMASS_Ecoli_core_w_GAM"]

# optimizer = OSQP.Optimizer 
# atts = Dict("eps_abs" => 5e-4,"eps_rel" => 5e-4, "max_iter" => 100_000, "verbose"=>false) # not a good linear solver :/
# sol = pfba(model, biomass, optimizer; solver_attributes=atts)
# sol["BIOMASS_Ecoli_core_w_GAM"]
