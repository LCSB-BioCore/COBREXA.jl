"""
    map_fluxes(v, model::CobraModel)

Map fluxes from an optimization problem (`v`) to rxns in a model.
`v` can be a JuMP object (fluxes) or an array of Float64 fluxes.
Assumes they are in order of `model.reactions`, which they should be since the optimization problem is constructed from the model.
"""
function map_fluxes(v::Array{Float64,1}, model::CobraModel)
    rxndict = Dict{String,Float64}()
    for i in eachindex(model.reactions)
        rxndict[model.reactions[i].id] = v[i]
    end
    return rxndict
end

function map_fluxes(v::Array{VariableRef,1}, model::CobraModel)
    rxndict = Dict{String,Float64}()
    for i in eachindex(model.reactions)
        rxndict[model.reactions[i].id] = value(v[i])
    end
    return rxndict
end

"""
    atom_exchange(flux_dict::Dict{String, Float64}, model::CobraModel)

Return a dictionary mapping the flux of atoms across the boundary of the model given `flux_dict` of reactions in `model`.
Here `flux_dict` is a mapping of reaction `id`s to fluxes, e.g. from FBA.
"""
function atom_exchange(flux_dict::Dict{String,Float64}, model::CobraModel)
    atom_flux = Dict{String,Float64}()
    for (rxnid, flux) in flux_dict
        if startswith(rxnid, "EX_") || startswith(rxnid, "DM_") # exchange, demand reaction
            for (met, stoich) in findfirst(model.reactions, rxnid).metabolites
                adict = get_atoms(met)
                for (atom, stoich) in adict
                    atom_flux[atom] = get(atom_flux, atom, 0.0) + flux * stoich
                end
            end
        end
    end
    return atom_flux
end

"""
    get_exchanges(rxndict::Dict{String, Float64}; top_n=8, ignorebound=1000.0, verbose=true)

Display the top_n producing and consuming exchange fluxes.
Set top_n to a large number to get all the consuming/producing fluxes.
Ignores infinite (problem upper/lower bound) fluxes (set with ignorebound).
When `verbose` is false, the output is not printed out.
Return these reactions in two dictionaries: `consuming`, `producing`
"""
function exchange_reactions(
    rxndict::Dict{String,Float64};
    top_n = 8,
    ignorebound = 1000.0,
    verbose = true,
)
    fluxes = Float64[]
    rxns = String[]
    for (k, v) in rxndict
        if startswith(k, "EX_") && abs(v) < ignorebound
            push!(rxns, k)
            push!(fluxes, v)
        end
    end
    inds_prod = sortperm(fluxes, rev = true)
    inds_cons = sortperm(fluxes)

    consuming = Dict{String,Float64}()
    producing = Dict{String,Float64}()
    verbose && println("Consuming fluxes:")
    for i = 1:min(top_n, length(rxndict))
        if rxndict[rxns[inds_cons[i]]] < -eps()
            verbose && println(
                rxns[inds_cons[i]],
                " = ",
                round(rxndict[rxns[inds_cons[i]]], digits = 4),
            )
            consuming[rxns[inds_cons[i]]] = rxndict[rxns[inds_cons[i]]]
        else
            continue
        end
    end

    verbose && println("Producing fluxes:")
    for i = 1:min(top_n, length(rxndict))
        if rxndict[rxns[inds_prod[i]]] > eps()
            verbose && println(
                rxns[inds_prod[i]],
                " = ",
                round(rxndict[rxns[inds_prod[i]]], digits = 4),
            )
            producing[rxns[inds_prod[i]]] = rxndict[rxns[inds_prod[i]]]
        else
            continue
        end
    end
    return consuming, producing
end

"""
    metabolite_fluxes(fluxdict::Dict{String, Float64}, model::CobraModel)

Return two dictionaries of metabolite `id`s mapped to reactions that consume or produce them given the flux distribution supplied in `fluxdict`.
"""
function metabolite_fluxes(fluxdict::Dict{String,Float64}, model::CobraModel)
    S = Array(stoichiometry(model))
    met_flux = Dict{String,Float64}()
    rxnids = reactions(model)
    metids = metabolites(model)

    producing = Dict{String,Dict{String,Float64}}()
    consuming = Dict{String,Dict{String,Float64}}()
    for (row, metid) in enumerate(metids)
        for (col, rxnid) in enumerate(rxnids)
            mf = fluxdict[rxnid] * S[row, col]
            # ignore zero flux
            if mf < -eps() # consuming rxn
                if haskey(consuming, metid)
                    consuming[metid][rxnid] = mf
                else
                    consuming[metid] = Dict(rxnid => mf)
                end
            elseif mf > eps()
                if haskey(producing, metid)
                    producing[metid][rxnid] = mf
                else
                    producing[metid] = Dict(rxnid => mf)
                end
            end
        end
    end
    return consuming, producing
end
