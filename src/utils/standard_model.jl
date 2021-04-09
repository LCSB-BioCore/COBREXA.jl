
"""
    getindex(model::StandardModel, rxn::Reaction)

Get the index of `rxn` in `model`, based on reaction `id`.
Return -1 if not found.

Typical usage: ind = model[rxn]
"""
function Base.getindex(model::StandardModel, rxn::Reaction)
    return model.reactions[rxn]
end

"""
    getindex(model::StandardModel, met::Metabolite)

Get the index of `met` in `model`, based on metabolite `id`.
Return -1 if not found.

Typical usage: ind = model[met]
"""
function Base.getindex(model::StandardModel, met::Metabolite)
    return model.metabolites[met]
end

"""
    getindex(model::StandardModel, gene::Gene)

Get the index of `gene` in `model`, based on gene `id`.
Return -1 if not found.

Typical usage: ind = model[gene]
"""
function Base.getindex(model::StandardModel, gene::Gene)
    return model.genes[gene]
end

"""
    atom_exchange(flux_dict::Dict{String, Float64}, model::StandardModel)

Return a dictionary mapping the flux of atoms across the boundary of the model given `flux_dict` of reactions in `model`.
Here `flux_dict` is a mapping of reaction `id`s to fluxes, e.g. from FBA.
"""
function atom_exchange(flux_dict::Dict{String,Float64}, model::StandardModel)
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
    get_exchanges(rxndict::Dict{String, Float64}; top_n=8, ignorebound=_constants.default_reaction_bound, verbose=true)

Display the top_n producing and consuming exchange fluxes.
Set top_n to a large number to get all the consuming/producing fluxes.
Ignores infinite (problem upper/lower bound) fluxes (set with ignorebound).
When `verbose` is false, the output is not printed out.
Return these reactions in two dictionaries: `consuming`, `producing`
"""
function exchange_reactions(
    rxndict::Dict{String,Float64};
    top_n = 8,
    ignorebound = _constants.default_reaction_bound,
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
    metabolite_fluxes(fluxdict::Dict{String, Float64}, model::StandardModel)

Return two dictionaries of metabolite `id`s mapped to reactions that consume or 
produce them given the flux distribution supplied in `fluxdict`.
"""
function metabolite_fluxes(fluxdict::Dict{String,Float64}, model::StandardModel)
    S = stoichiometry(model)
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

"""
    set_bound(index, optimization_model;
        ub=_constants.default_reaction_rate,
        lb=-_constants.default_reaction_rate)

Helper function to set the bounds of variables.
The JuMP `set_normalized_rhs` function is a little confusing, 
so this function simplifies setting constraints. In short, JuMP
uses a normalized right hand side representation of constraints, 
which means that lower bounds have their sign flipped. This function
does this for you, so you don't have to remember to do this whenever you
change the constraints. 

Just supply the constraint `index` and the JuMP model (`opt_model`) that 
will be solved, and the variable's bounds will be set to `ub` and `lb`.
"""
function set_bound(
    vind,
    opt_model;
    ub = _constants.default_reaction_rate,
    lb = -_constants.default_reaction_rate,
)
    if lb <= 0
        set_normalized_rhs(opt_model[:lbs][vind], abs(lb))
    else
        set_normalized_rhs(opt_model[:lbs][vind], -abs(lb))
    end
    set_normalized_rhs(opt_model[:ubs][vind], ub)
end

"""
    get_bound_vectors(opt_model)

Returns vectors of the lower and upper bounds of `opt_model` constraints, where 
`opt_model` is a JuMP model constructed by e.g. `make_optimization_problem` or
`flux_balance_analysis`.

See also: [`make_optimization_problem`](@ref), [`flux_balance_analysis`](`ref`)
"""
function get_bound_vectors(opt_model)
    lbconref = opt_model[:lbs]
    ubconref = opt_model[:ubs]
    lbs = zeros(length(lbconref))
    for i in eachindex(lbs)
        lbval = normalized_rhs(lbconref[i])
        if lbval > 0
            lbs[i] = -abs(lbval)
        else
            lbs[i] = abs(lbval)
        end
    end
    ubs = [normalized_rhs(ubconref[i]) for i in eachindex(ubconref)]

    return lbs, ubs
end
