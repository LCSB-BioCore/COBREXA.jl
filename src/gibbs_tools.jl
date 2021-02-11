# mutable struct Gibbs
#     rxn :: Reaction
#     ΔG :: Measurement{Float64}
# end
# mutable struct ReactionGibbs
    
# end

# function Base.show(io::IO, rgs::ReactionGibbs)
#     println(io, "Gibbs free energy of reaction set of length: $(length(rgs.rxns))")
# end

# function Base.getindex(rgs::ReactionGibbs, rxn::Reaction)
#     for 
#     println(io, "Reaction ID: ")
#     println(io, "ΔG = ", rgs)
# end

# function mapGibbs(rxns::Array{Reaction, 1}; dgtype="zero", ph=7.0, ionic_str="100 mM") 
#     gibbs = Gibbs[]
#     for rxn in rxns
#         # println(rxn.id)
#         f = buildrxnstring(rxn)
#         try 
#             if dgtype == "phys"
#                 ΔG = Gibbs(getdG0(f, ph, ionic_str))
#             elseif dgtype == "prime" 
#                 ΔG = Gibbs(getdGprime(f, ph, ionic_str))
#             else # "zero"
#                 ΔG = Gibbs(getdGphys(f, ph, ionic_str))
#             end

#             push!(gibbs, ΔG)
#         catch
#             @warn "Error on reaction: $(rxn.id)"
#         end
#     end
#     return ReactionGibbs(rxns, gibbs)
# end

# function getdG0(formula::String, ph=7.0, ionic_strength="100 mM")
#     isbal, v, err = py"pygetdg0"(formula, ph, ionic_strength)
#     !isbal && @warn "Reaction not balanced: $formula"
#     return v ± err
# end

# function getdGprime(formula::String, ph=7.0, ionic_strength="100 mM")
#     isbal, v, err = py"pygetdgprime"(formula, ph, ionic_strength)
#     !isbal && @warn "Reaction not balanced: $formula"
#     return v ± err
# end

# function getdGphys(formula::String, ph=7.0, ionic_strength="100 mM")
#     isbal, v, err = py"pygetdgprimephys"(formula, ph, ionic_strength)
#     !isbal && @warn "Reaction not balanced: $formula"
#     return v ± err
# end

# function blackbox(fluxres, gibbs, biomassrxn, biomassdg)
#     total_ΔG = 0.0 ± 0.0
#     for rxn in fluxres.rxns
#         if startswith(rxn.id, "EX_")
#             gind = gibbs.rxns[rxn]
#             fluxind = fluxres.rxns[rxn]
#             total_ΔG += fluxres.fluxes[fluxind] * gibbs.ΔGs[gind].ΔG
#         end
#     end
#     total_ΔG += fluxres.fluxes[fluxres.rxns[biomassrxn]]*biomassdg
#     return total_ΔG
# end

# """
# buildrxnstring(model, rxnid)

# Get rxn in string format for Equilibrator.
# """
# function buildrxnstring(rxn::Reaction)
#     pos_s = []
#     neg_s = []
#     for (met, coeff) in rxn.metabolites 
#         if coeff > 0.0
#             push!(pos_s, "$(coeff) bigg.metabolite:$(met.id[1:end-2])")
#         else
#             push!(neg_s, "$(abs(coeff)) bigg.metabolite:$(met.id[1:end-2])")    
#         end
#     end
#     return join(neg_s, " + ")*" = "*join(pos_s, " + ") # keep order for ease of use later
# end