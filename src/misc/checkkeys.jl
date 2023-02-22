"""
Throw a DomainError if `model.field` has any keys in the ID field of
xs::Union{Reaction, Metabolite, Gene}.
"""
function check_arg_keys_exists(model, field, xs)
    d = getfield(model, field)
    ids = [x.id for x in xs if haskey(d, x.id)]
    isempty(ids) ||
        throw(DomainError("Duplicated $field IDs already present in model: $ids"))
    nothing
end

"""
Throw a DomainError if `model.field` does not have any keys in `xs`.
"""
function check_arg_keys_missing(model, field, xs::Vector{String})
    d = getfield(model, field)
    ids = [x for x in xs if !haskey(d, x)]
    isempty(ids) || throw(DomainError(ids, " $field IDs not found in model."))
    nothing
end

"""
Throw a DomainError if rxn_ids have isozymes associated with them.
"""
function check_has_isozymes(model, rxn_ids)
    ids = [rid for rid in rxn_ids if !isnothing(model.reactions[rid].gene_associations)]
    isempty(ids) || throw(DomainError(ids, " reactions already have isozymes."))
    nothing
end

"""
Throw a DomainError if the `biomass_rxn_id` is not in `model_reactions`.
"""
function check_has_biomass_rxn_id(model_reactions, biomass_rxn_id)
    haskey(model_reactions, biomass_rxn_id) ||
        throw(DomainError(biomass_rxn_id, " not found in model."))
    nothing
end

"""
Throw a DomainError if the `biomass_rxn_id` in `model_reactions` has any
isozymes assigned to it.
"""
function check_biomass_rxn_has_isozymes(model_reactions, biomass_rxn_id)
    isnothing(model_reactions[biomass_rxn_id].gene_associations) ||
        throw(DomainError(biomass_rxn_id, " already has isozymes associated to it."))
    nothing
end

"""
Throw a DomainError if `virtualribosome_id` is already in the `model_genes`.
"""
function check_has_virtualribosome(model_genes, virtualribosome_id)
    haskey(model_genes, virtualribosome_id) &&
        throw(DomainError(virtualribosome_id, " already found in model."))
    nothing
end

"""
Throw a DomainError if `biomass_rxn_id` in `model_reactions` already has a
`biomass_metabolite_id`.
"""
function check_has_biomass_rxn_biomas_metabolite(
    model_reactions,
    biomass_rxn_id,
    biomass_metabolite_id,
)
    haskey(model_reactions[biomass_rxn_id], biomass_metabolite_id) ||
        throw(DomainError(biomass_metabolite_id, " not found in $biomass_rxn_id."))
    nothing
end

"""
Throw a DomainError if some reaction ids `rids` are not found in the
environmental_links of `cm`. Return the indices of `rids` in the environmental
linkage vector.
"""
function check_environmental_ids(cm, rids)
    env_rids = [envlink.reaction_id for envlink in cm.environmental_links]
    idxs = indexin(rids, env_rids)
    any(isnothing.(idxs)) || begin 
        missing_idxs = findall(isnothing, idxs)
        throw(DomainError(rids[missing_idxs], " exchange reaction IDs not found in environmental links."))
    end
    idxs
end
