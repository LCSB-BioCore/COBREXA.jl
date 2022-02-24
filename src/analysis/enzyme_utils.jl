"""
    _bounds(model::StandardModel, rid::String)

Return lower and upper bounds for `rid` in `model`.
"""
function _bounds(model::StandardModel, rid::String)
    #TODO generalize this to other model types
    model.reactions[rid].lb, model.reactions[rid].ub
end

"""
    _is_reversible(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is reversible.
"""
function _is_reversible(model::StandardModel, rid::String)
    lb, ub = _bounds(model, rid)
    lb < 0 && ub > 0
end

"""
    _is_forward_only(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is forward only.
"""
function _is_forward_only(model::StandardModel, rid::String)
    lb, ub = _bounds(model, rid)
    lb >= 0 && ub > 0
end

"""
    _is_backward_only(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is backward only.
"""
function _is_backward_only(model::StandardModel, rid::String)
    lb, ub = _bounds(model, rid)
    lb < 0 && ub <= 0
end

"""
    _is_unidirectional(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is unidirectional.
"""
function _is_unidirectional(model::StandardModel, rid::String)
    _is_forward_only(model, rid) || _is_backward_only(model, rid)
end

"""
    _is_blocked(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is blocked.
"""
function _is_blocked(model::StandardModel, rid::String)
    lb, ub = _bounds(model, rid)
    lb == ub == 0
end

"""
    _has_isozymes(model::StandardModel, rid::String)

Check if reaction `rid` in `model` is catalyzed by multiple enzymes, 
i.e. it has isozymes according to the gene reaction rules.
"""
function _has_isozymes(model::StandardModel, rid::String)
    length(reaction_gene_association(model, rid)) > 1
end

"""
    _has_grr(model::StandardModel, rid::String)

Check if reaction `rid` in `model` has a gene reaction rule entry.
"""
function _has_grr(model::StandardModel, rid::String)
    #TODO simplify this once COBREXA enforces universal rules for GRR representation
    !isnothing(reaction_gene_association(model, rid)) &&
        reaction_gene_association(model, rid) != [[]] &&
        !isempty(first(reaction_gene_association(model, rid)))
end

"""
    _get_proteins_with_kcats(model::StandardModel, reaction_kcats)

Return all protein (gene ids) that have a kcat from `model` based on `reaction_kcats`, 
which is a dictionary mapping reaction ids to the kcats of each isozyme. Assume that if 
a reaction has a kcat then each isozyme has a kcat.
"""
function _get_proteins_with_kcats(model::StandardModel, reaction_kcats)
    unique(
        vcat(
            vcat(
                [
                    reaction_gene_association(model, rid) for
                    rid in reactions(model) if haskey(reaction_kcats, rid)
                ]...,
            )...,
        ),
    )
end

"""
    _build_irreversible_stoichiometric_matrix(model::StandardModel)

Return the stoichiometric matrix. All reactions are forward only i.e. only 
positive fluxes are allowed. Include arm reactions.
"""
function _build_irreversible_stoichiometric_matrix(model::StandardModel)
    # components used to build stoichiometric matrix
    S_components = ( #TODO add size hints if possible
        row_idxs = Vector{Int}(),
        col_idxs = Vector{Int}(),
        coeffs = Vector{Float64}(),
        lbs = Vector{Float64}(),
        ubs = Vector{Float64}(),
    )

    # establish the ordering in a named tuple
    idxs = ( #: pseudo metabolites and reactions are added to model
        met_idxs = Dict{String,Int}(),
        rxn_idxs = Dict{String,Int}(),
        max_rxn_idx = [1], #TODO maybe fix, this is a dodgy way of adding a counter to a named tuple
        max_met_idx = [1], #TODO maybe fix, this is a dodgy way of adding a counter to a named tuple
        pseudo_met_idx = [1], #TODO maybe fix, this is a dodgy way of adding a counter to a named tuple
    )
    #TODO for the counter thing, basically I wanted e.g. max_rxn_idx = 1 and then update it, 
    #TODO but named tuples are immutable... :(

    # fill the matrix entries
    #: blocked treated as reversible because unclear what direction the reaction would go
    for rid in reactions(model)
        if _has_grr(model, rid) && _has_isozymes(model, rid)
            if _is_unidirectional(model, rid)
                dir = _is_forward_only(model, rid) ? "§FOR" : "§REV"
                _add_isozyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)
            elseif _is_reversible(model, rid) || _is_blocked(model, rid)
                _add_isozyme_to_irrev_stoich_mat(model, rid, idxs, S_components, "§FOR")
                _add_isozyme_to_irrev_stoich_mat(model, rid, idxs, S_components, "§REV")
            else
                @warn "Unhandled bound type for $rid"
            end
        else # no grr or single enzyme only
            if _is_unidirectional(model, rid)
                dir = _is_forward_only(model, rid) ? "§FOR" : "§REV"
                _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)
            elseif _is_reversible(model, rid) || _is_blocked(model, rid)
                _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, "§FOR")
                _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, "§REV")
            else
                @warn "Unhandled bound type for $rid"
            end
        end
    end

    S = sparse(
        S_components.row_idxs,
        S_components.col_idxs,
        S_components.coeffs,
        length(idxs.met_idxs),
        length(idxs.rxn_idxs),
    )

    return S, S_components.lbs, S_components.ubs, idxs.rxn_idxs, idxs.met_idxs
end


"""
    _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)

Add entries to the components that will be used to build the stoichiometric matrix.
Simple variant that does not deal with isozymes and arm reactions.
"""
function _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)
    idxs.rxn_idxs[rid*dir] = idxs.max_rxn_idx[1]
    idxs.max_rxn_idx[1] += 1
    fix_sign = dir == "§FOR" ? 1 : -1 # change direction of reaction
    for (mid, coeff) in reaction_stoichiometry(model, rid)
        if !haskey(idxs.met_idxs, mid)
            idxs.met_idxs[mid] = idxs.max_met_idx[1]
            idxs.max_met_idx[1] += 1
        end
        push!(S_components.row_idxs, idxs.met_idxs[mid])
        push!(S_components.col_idxs, idxs.rxn_idxs[rid*dir])
        push!(S_components.coeffs, fix_sign * coeff)
    end
    lb, ub = abs.(_bounds(model, rid)) # assumes lb < ub
    if dir == "§FOR"
        _is_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, lb)
        push!(S_components.ubs, ub)
    else
        _is_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, ub)
        push!(S_components.ubs, lb)
    end
end

"""
    _add_isozyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)

Add entries to the components that will be used to build the stoichiometric matrix.
Complex variant that deals with isozymes and arm reactions.
"""
function _add_isozyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)
    # add pseudo metabolite
    pm = "§PM$(idxs.pseudo_met_idx[1])"
    idxs.pseudo_met_idx[1] += 1
    idxs.met_idxs[pm] = idxs.max_met_idx[1]
    idxs.max_met_idx[1] += 1
    # find half reactions to get arm reaction
    lhs = []
    rhs = []
    for (mid, coeff) in reaction_stoichiometry(model, rid)
        if !haskey(idxs.met_idxs, mid)
            idxs.met_idxs[mid] = idxs.max_met_idx[1]
            idxs.max_met_idx[1] += 1
        end
        if coeff <= 0
            push!(lhs, (mid, coeff))
        else
            push!(rhs, (mid, coeff))
        end
    end
    product_half_reaction = dir == "§FOR" ? rhs : lhs
    reagent_half_reaction = dir == "§FOR" ? lhs : rhs
    # add arm reaction
    fix_sign = dir == "§FOR" ? 1 : -1 # change direction of reaction
    pr = rid * "§ARM" * dir
    idxs.rxn_idxs[pr] = idxs.max_rxn_idx[1] #! this needs to get added first because of blocked possibility
    idxs.max_rxn_idx[1] += 1
    push!(S_components.row_idxs, idxs.met_idxs[pm])
    push!(S_components.col_idxs, idxs.rxn_idxs[pr])
    push!(S_components.coeffs, 1)
    for (mid, coeff) in reagent_half_reaction
        push!(S_components.row_idxs, idxs.met_idxs[mid])
        push!(S_components.col_idxs, idxs.rxn_idxs[pr])
        push!(S_components.coeffs, fix_sign * coeff)
    end
    # add bounds for ARM reaction that corresponds to original model's bounds
    lb, ub = abs.(_bounds(model, rid)) # assumes lb < ub
    if dir == "§FOR"
        _is_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, lb)
        push!(S_components.ubs, ub)
    else
        _is_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, ub)
        push!(S_components.ubs, lb)
    end
    # add isozyme reactions
    for (i, _) in enumerate(reaction_gene_association(model, rid))
        iso_rid = rid * "§ISO$i" * dir
        idxs.rxn_idxs[iso_rid] = idxs.max_rxn_idx[1]
        idxs.max_rxn_idx[1] += 1
        push!(S_components.row_idxs, idxs.met_idxs[pm])
        push!(S_components.col_idxs, idxs.rxn_idxs[iso_rid])
        push!(S_components.coeffs, -1)
        for (mid, coeff) in product_half_reaction
            push!(S_components.row_idxs, idxs.met_idxs[mid])
            push!(S_components.col_idxs, idxs.rxn_idxs[iso_rid])
            push!(S_components.coeffs, fix_sign * coeff)
        end
        # add bounds
        push!(S_components.lbs, 0)
        if _is_blocked(model, rid)
            push!(S_components.ubs, 0)
        else
            push!(S_components.ubs, 10_000) # arbitrary upper bound
        end
    end
end

"""
    _add_enzyme_variable(model, iso_num, rid, original_rid, protein_stoichiometry, reaction_kcats, E_components, col_idx, protein_ids)

Helper function to add an column into the enzyme stoichiometric matrix.
"""
function _add_enzyme_variable(
    model,
    iso_num,
    rid,
    original_rid,
    protein_stoichiometry,
    reaction_kcats,
    E_components,
    col_idx,
    protein_ids,
)
    grr = reaction_gene_association(model, original_rid)[iso_num]
    pstoich = protein_stoichiometry[original_rid][iso_num]
    kcat =
        contains(rid, "§FOR") ? reaction_kcats[original_rid][iso_num][1] :
        reaction_kcats[original_rid][iso_num][2]
    for (idx, pid) in enumerate(grr)
        push!(E_components.row_idxs, first(indexin([pid], protein_ids)))
        push!(E_components.col_idxs, col_idx)
        push!(E_components.coeffs, -pstoich[idx] / kcat)
    end
end

"""
    _order_id_to_idx_dict(id_to_idx_dict)

Return the keys of `id_to_idx_dict` sorted by the values, which
are taken to be the indices.
"""
function _order_id_to_idx_dict(dmap)
    ks = collect(keys(dmap))
    vs = collect(values(dmap))
    return ks[sortperm(vs)]
end

"""
    _map_irrev_to_rev_ids(reaction_map, protein_ids, solution)

Return dictionaries of reaction ids mapped to fluxes, 
and protein ids mapped to concentrations using `reaction_map` to 
determine the ids of fluxes and `protein_ids` for the gene ids.
The solution in `solution` is used to fill the dictionaries.
"""
function _map_irrev_to_rev_ids(reaction_map, solution; protein_ids = [])
    reaction_flux = Dict{String,Float64}()
    for (k, i) in reaction_map
        contains(k, "§ISO") && continue # §ISO§FOR and §ISO§REV need to be ignored
        rid = split(k, "§")[1]
        v = contains(k, "§FOR") ? solution[i] : -solution[i]
        reaction_flux[rid] = get(reaction_flux, rid, 0) + v
    end

    if isempty(protein_ids)
        return reaction_flux
    else
        n_reactions = length(reaction_map)
        protein_flux = Dict{String,Float64}()
        for (i, pid) in enumerate(protein_ids)
            protein_flux[pid] = solution[n_reactions+i]
        end
        return reaction_flux, protein_flux
    end
end

"""
    remove_slow_isozymes!(
        model::StandardModel;
        reaction_kcats = Dict(),
        protein_stoichiometry = Dict(),
        protein_masses = Dict(),
    )

Remove all but the fastest isozyme from each reaction in `model`.
Use the largest kcat (for, rev) for these calculations. Modifies all 
the arguments in place.
"""
function remove_slow_isozymes!(
    model::StandardModel;
    reaction_kcats = Dict(),
    protein_stoichiometry = Dict(),
    protein_masses = Dict(),
)
    for rid in reactions(model)
        if _has_grr(model, rid) && haskey(reaction_kcats, rid)
            kcat_effs = Float64[]
            grrs = reaction_gene_association(model, rid)
            for (i, grr) in enumerate(grrs)
                push!(
                    kcat_effs,
                    dot(
                        protein_stoichiometry[rid][i],
                        [protein_masses[gid] for gid in grr],
                    ) / maximum(reaction_kcats[rid][i]),
                )
            end
            idx = argmin(kcat_effs)

            model.reactions[rid].grr = [grrs[idx]]
            reaction_kcats[rid] = [reaction_kcats[rid][idx]]
            protein_stoichiometry[rid] = [protein_stoichiometry[rid][idx]]
        end
    end

    curated_gids = String[]
    for rid in reactions(model)
        if _has_grr(model, rid)
            for grr in reaction_gene_association(model, rid)
                append!(curated_gids, grr)
            end
        end
    end
    rm_gids = setdiff(genes(model), curated_gids)
    delete!(model.genes, rm_gids) # remove genes that were deleted

    return nothing
end
