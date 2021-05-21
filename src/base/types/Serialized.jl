
"""
    mutable struct Serialized{M <: MetabolicModel}
        m::Maybe{M}
        filename::String
    end

A meta-model that represents a model that is serialized on the disk. The
internal model will be loaded on-demand by using any accessor, or by calling
[`precache!`](@ref) directly.
"""
mutable struct Serialized{M} <: MetabolicModel where {M<:MetabolicModel}
    m::Maybe{M}
    filename::String
end

function _on_precached(m::Serialized, f)
    precache!(m)
    f(m.m)
end

reactions(m::Serialized) = _on_precached(m, reactions)
n_reactions(m::Serialized) = _on_precached(m, n_reactions)
metabolites(m::Serialized) = _on_precached(m, metabolites)
n_metabolites(m::Serialized) = _on_precached(m, n_metabolites)
stoichiometry(m::Serialized) = _on_precached(m, stoichiometry)
bounds(m::Serialized) = _on_precached(m, bounds)
balance(m::Serialized) = _on_precached(m, balance)
objective(m::Serialized) = _on_precached(m, objective)
coupling(m::Serialized) = _on_precached(m, coupling)
n_coupling_constraints(m::Serialized) = _on_precached(m, n_coupling_constraints)
coupling_bounds(m::Serialized) = _on_precached(m, coupling_bounds)
genes(m::Serialized) = _on_precached(m, genes)
n_genes(m::Serialized) = _on_precached(m, n_genes)
metabolite_formula(m::Serialized) = _on_precached(m, metabolite_formula)
metabolite_charge(m::Serialized) = _on_precached(m, metabolite_charge)
reaction_annotations(m::Serialized) = _on_precached(m, reaction_annotations)
metabolite_annotations(m::Serialized) = _on_precached(m, metabolite_annotations)
gene_annotations(m::Serialized) = _on_precached(m, gene_annotations)
reaction_nodes(m::Serialized) = _on_precached(m, reaction_nodes)
metabolite_nodes(m::Serialized) = _on_precached(m, metabolite_nodes)
gene_notes(m::Serialized) = _on_precached(m, gene_notes)
metabolite_compartment(m::Serialized) = _on_precached(m, metabolite_compartment)
reaction_subsystem(m::Serialized) = _on_precached(m, reaction_subsystem)

"""
    precache!(model::Serialized{MetabolicModel})::Nothing

Load the `Serialized` model from disk in case it's not alreadly loaded.
"""
function precache!(model::Serialized)::Nothing
    if isnothing(model.m)
        model.m = deserialize(model.filename)
    end
    nothing
end
