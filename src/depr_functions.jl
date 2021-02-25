# """
# ReactionFluxes holds references to the objective, reactions and fluxes from analysis (e.g. FBA)
# """
# mutable struct ReactionFluxes
#     objective_id :: String 
#     objective :: Float64
#     rxns :: Array{Reaction, 1}
#     fluxes :: Array{Float64, 1}
# end

# """
# getindex(reactionfluxes, rxn)

# Return the index of rxn in reactionfluxes and -1 if it is not found.
# Note, this is slightly different from the normal getindex function.
# """
# function Base.getindex(rfs::ReactionFluxes, rxn::Reaction)
#     for i in eachindex(rfs.rxns)
#         if rxn.id == rfs.rxns[i].id
#             return i
#         end
#     end
#     return -1
# end

# """
# Pretty printing of ReactionFluxes objects.
# """
# function Base.show(io::IO, rfs::ReactionFluxes)
#     println(io, "Optimum for $(rfs.objective_id) = ", round(rfs.objective, digits=4))
    
#     println()
#     inds = sortperm(rfs.fluxes) # consuming fluxes  
#     println("Consuming fluxes:")
#     counter = 0
#     for i in inds
#         if startswith(rfs.rxns[i].id, "EX_")
#             println(io, rfs.rxns[i].name, " = ", round(rfs.fluxes[i], digits=4), " mmol/gDW/h")
#             counter += 1
#         end
#         if counter > 8 # only display top 10
#             break
#         end
#     end

#     println()
#     inds = sortperm(rfs.fluxes, rev=true) # consuming fluxes  
#     println("Producing fluxes:")
#     counter = 0
#     for i in inds
#         if startswith(rfs.rxns[i].id, "EX_")
#             println(io, rfs.rxns[i].name, " = ", round(rfs.fluxes[i], digits=4), " mmol/gDW/h")
#             counter += 1
#         end
#         if counter > 8 # only display top 10
#             break
#         end
#     end
# end

# function atom_exchange(rfs::ReactionFluxes)
#     # find exchange reactions
#     ex_inds = [i for i in eachindex(rfs.rxns) if startswith(rfs.rxns[i].id, "EX_")]
    
#     atom_balance = Dict{String, Float64}()
#     for ex_ind in ex_inds
#         for (met, w) in rfs.rxns[ex_ind].metabolites
#             for (atom, stoich) in get_atoms(met)
#                 atom_balance[atom] = get(atom_balance, atom, 0.0) + stoich*w*value(rfs.fluxes[ex_ind])
#             end
#         end
#     end
#     return atom_balance
# end

# function map_gibbs_external(fluxres::ReactionFluxes, gibbs)
#     total_ΔG = 0.0 ± 0.0
#     missing_flux = 0.0
#     for (i, rxn) in enumerate(fluxres.rxns)
#         if startswith(rxn.id, "EX_")
#             if gibbs[rxn.id] ≈ 0.0
#                 missing_flux += abs(fluxres.fluxes[i])
#             end 
#             total_ΔG -= fluxres.fluxes[i] * gibbs[rxn.id] # negative here because formation is not MET -> as used here, but the -> MET 
#         end
#     end
#     return total_ΔG, missing_flux/sum(abs, fluxres.fluxes) # units J/gDW/h
# end

# function map_gibbs_internal(fluxres::ReactionFluxes, gibbs, biomassid="BIOMASS")
#     total_ΔG = 0.0 ± 0.0
#     missing_flux = 0.0
#     found_flux = 0.0
#     for (i, rxn) in enumerate(fluxres.rxns)
#         if !startswith(rxn.id, "EX_") && !contains(rxn.id, biomassid) # ignore exchange reactions and biomass eqn
#             if gibbs[rxn.id] ≈ 0.0
#                 missing_flux += abs(fluxres.fluxes[i])
#             else
#                 found_flux += abs(fluxres.fluxes[i])
#             end 
#             total_ΔG += fluxres.fluxes[i] * gibbs[rxn.id] # add because this is not formation but rather just adding equations (the flux direction sign compensates) 
#         end
#     end
#     return total_ΔG, missing_flux/(missing_flux+found_flux) # units J/gDW/h
# end
