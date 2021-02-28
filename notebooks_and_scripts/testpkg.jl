using CobraTools
using JuMP
using Gurobi
using Tulip
using OSQP
using GLPK
using Gadfly


# using Compose # export plot to file
# using Gadfly # export plot to browser
# using LightGraphs
# using GraphPlot

# g = Graph(4)
# add_edge!(g,1,2)
# add_edge!(g,1,3)
# add_edge!(g,2,4)
# add_edge!(g,1,4)
# gplot(g)

# draw(PNG("mygraph.png", 8cm, 8cm), gplot(g))


# optimizer = Tulip.Optimizer # quiet by default
# sol = fba(model, biomass, optimizer)
# sol["BIOMASS_Ecoli_core_w_GAM"]


# atts = Dict("eps_abs" => 5e-4,"eps_rel" => 5e-4, "max_iter" => 100_000, "verbose"=>false) # not a good linear solver :/
# sol = pfba(model, biomass, [Tulip.Optimizer, OSQP.Optimizer]; solver_attributes=Dict("opt1" => Dict{Any, Any}(), "opt2" => atts))
# sol["BIOMASS_Ecoli_core_w_GAM"]

# biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
# optimizer = OSQP.Optimizer 
# atts = Dict("eps_abs" => 5e-4,"eps_rel" => 5e-4, "max_iter" => 100_000, "verbose"=>false) # not a good linear solver :/
# sol = pfba(model, biomass, optimizer; solver_attributes=atts)
# sol["PGM"]
