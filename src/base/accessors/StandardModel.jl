
# MetabolicModel interface follows
"""
$(TYPEDSIGNATURES)

Return a vector of reaction id strings contained in `model`.
The order of reaction ids returned here matches the order used to construct the
stoichiometric matrix.
"""
reactions(model::StandardModel)::StringVecType = collect(keys(model.reactions))

"""
$(TYPEDSIGNATURES)

Return the number of reactions contained in `model`.
"""
n_reactions(model::StandardModel)::Int = length(model.reactions)


"""
$(TYPEDSIGNATURES)

Return a vector of metabolite id strings contained in `model`.
The order of metabolite strings returned here matches the order used to construct
the stoichiometric matrix.
"""
metabolites(model::StandardModel)::StringVecType = collect(keys(model.metabolites))

"""
$(TYPEDSIGNATURES)

Return the number of metabolites in `model`.
"""
n_metabolites(model::StandardModel)::Int = length(model.metabolites)

"""
$(TYPEDSIGNATURES)

Return a vector of gene id strings in `model`.
"""
genes(model::StandardModel)::StringVecType = collect(keys(model.genes))

"""
$(TYPEDSIGNATURES)

Return the number of genes in `model`.
"""
n_genes(model::StandardModel)::Int = length(model.genes)

"""
$(TYPEDSIGNATURES)

Return the stoichiometric matrix associated with `model` in sparse format.
"""
function stoichiometry(model::StandardModel)::SparseMat
    n_entries = 0
    for (_, r) in model.reactions
        for _ in r.metabolites
            n_entries += 1
        end
    end

    MI = Vector{Int}()
    RI = Vector{Int}()
    SV = Vector{Float64}()
    sizehint!(MI, n_entries)
    sizehint!(RI, n_entries)
    sizehint!(SV, n_entries)

    # establish the ordering
    rxns = reactions(model)
    met_idx = Dict(mid => i for (i, mid) in enumerate(metabolites(model)))

    # fill the matrix entries
    for (ridx, rid) in enumerate(rxns)
        for (mid, coeff) in model.reactions[rid].metabolites
            haskey(met_idx, mid) || throw(
                DomainError(
                    mid,
                    "Metabolite $(mid) not found in model but occurs in stoichiometry of $(rid)",
                ),
            )
            push!(MI, met_idx[mid])
            push!(RI, ridx)
            push!(SV, coeff)
        end
    end
    return SparseArrays.sparse(MI, RI, SV, n_metabolites(model), n_reactions(model))
end

"""
$(TYPEDSIGNATURES)

Return the lower bounds for all reactions in `model` in sparse format.
"""
lower_bounds(model::StandardModel)::Vector{Float64} =
    sparse([model.reactions[rxn].lb for rxn in reactions(model)])

"""
$(TYPEDSIGNATURES)

Return the upper bounds for all reactions in `model` in sparse format.
Order matches that of the reaction ids returned in `reactions()`.
"""
upper_bounds(model::StandardModel)::Vector{Float64} =
    sparse([model.reactions[rxn].ub for rxn in reactions(model)])

"""
$(TYPEDSIGNATURES)

Return the lower and upper bounds, respectively, for reactions in `model`.
Order matches that of the reaction ids returned in `reactions()`.
"""
bounds(model::StandardModel)::Tuple{Vector{Float64},Vector{Float64}} =
    (lower_bounds(model), upper_bounds(model))

"""
$(TYPEDSIGNATURES)

Return the balance of the linear problem, i.e. b in Sv = 0 where S is the stoichiometric matrix
and v is the flux vector.
"""
balance(model::StandardModel)::SparseVec = spzeros(length(model.metabolites))

"""
$(TYPEDSIGNATURES)

Return sparse objective vector for `model`.
"""
objective(model::StandardModel)::SparseVec =
    sparse([model.reactions[rid].objective_coefficient for rid in keys(model.reactions)])

"""
$(TYPEDSIGNATURES)

Return the gene reaction rule in string format for reaction with `id` in `model`.
Return `nothing` if not available.
"""
reaction_gene_association(model::StandardModel, id::String)::Maybe{GeneAssociation} =
    _maybemap(identity, model.reactions[id].grr)

"""
$(TYPEDSIGNATURES)

Return the formula of reaction `id` in `model`.
Return `nothing` if not present.
"""
metabolite_formula(model::StandardModel, id::String)::Maybe{MetaboliteFormula} =
    _maybemap(_parse_formula, model.metabolites[id].formula)

"""
$(TYPEDSIGNATURES)

Return the charge associated with metabolite `id` in `model`.
Return nothing if not present.
"""
metabolite_charge(model::StandardModel, id::String)::Maybe{Int} =
    model.metabolites[id].charge

"""
$(TYPEDSIGNATURES)

Return compartment associated with metabolite `id` in `model`.
Return `nothing` if not present.
"""
metabolite_compartment(model::StandardModel, id::String)::Maybe{String} =
    model.metabolites[id].compartment

"""
$(TYPEDSIGNATURES)

Return the subsystem associated with reaction `id` in `model`.
Return `nothing` if not present.
"""
reaction_subsystem(model::StandardModel, id::String)::Maybe{String} =
    model.reactions[id].subsystem

"""
$(TYPEDSIGNATURES)

Return the notes associated with metabolite `id` in `model`.
Return an empty Dict if not present.
"""
metabolite_notes(model::StandardModel, id::String)::Maybe{Notes} =
    model.metabolites[id].notes

"""
$(TYPEDSIGNATURES)

Return the annotation associated with metabolite `id` in `model`.
Return an empty Dict if not present.
"""
metabolite_annotations(model::StandardModel, id::String)::Maybe{Annotations} =
    model.metabolites[id].annotations

"""
$(TYPEDSIGNATURES)

Return the notes associated with gene `id` in `model`.
Return an empty Dict if not present.
"""
gene_notes(model::StandardModel, id::String)::Maybe{Notes} = model.genes[id].notes

"""
$(TYPEDSIGNATURES)

Return the annotation associated with gene `id` in `model`.
Return an empty Dict if not present.
"""
gene_annotations(model::StandardModel, id::String)::Maybe{Annotations} =
    model.genes[id].annotations

"""
$(TYPEDSIGNATURES)

Return the notes associated with reaction `id` in `model`.
Return an empty Dict if not present.
"""
reaction_notes(model::StandardModel, id::String)::Maybe{Notes} = model.reactions[id].notes

"""
$(TYPEDSIGNATURES)

Return the annotation associated with reaction `id` in `model`.
Return an empty Dict if not present.
"""
reaction_annotations(model::StandardModel, id::String)::Maybe{Annotations} =
    model.reactions[id].annotations

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction with ID `rid`.
"""
reaction_stoichiometry(m::StandardModel, rid::String)::Dict{String,Float64} =
    m.reactions[rid].metabolites

"""
$(TYPEDSIGNATURES)

Return the name of reaction with ID `id`.
"""
reaction_name(m::StandardModel, rid::String) = m.reactions[rid].name

"""
$(TYPEDSIGNATURES)

Return the name of metabolite with ID `id`.
"""
metabolite_name(m::StandardModel, mid::String) = m.metabolites[mid].name

"""
$(TYPEDSIGNATURES)

Return the name of gene with ID `id`.
"""
gene_name(m::StandardModel, gid::String) = m.genes[gid].name

"""
$(TYPEDSIGNATURES)

Convert any `MetabolicModel` into a `StandardModel`.
Note, some data loss may occur since only the generic interface is used during
the conversion process.
"""
function Base.convert(::Type{StandardModel}, model::MetabolicModel)
    if typeof(model) == StandardModel
        return model
    end

    id = "" # TODO: add accessor to get model ID
    modelreactions = OrderedDict{String,Reaction}()
    modelmetabolites = OrderedDict{String,Metabolite}()
    modelgenes = OrderedDict{String,Gene}()

    gids = genes(model)
    metids = metabolites(model)
    rxnids = reactions(model)

    for gid in gids
        modelgenes[gid] = Gene(
            gid;
            name = gene_name(model, gid),
            notes = gene_notes(model, gid),
            annotations = gene_annotations(model, gid),
        )
    end

    for mid in metids
        modelmetabolites[mid] = Metabolite(
            mid;
            name = metabolite_name(model, mid),
            charge = metabolite_charge(model, mid),
            formula = _maybemap(_unparse_formula, metabolite_formula(model, mid)),
            compartment = metabolite_compartment(model, mid),
            notes = metabolite_notes(model, mid),
            annotations = metabolite_annotations(model, mid),
        )
    end

    S = stoichiometry(model)
    lbs, ubs = bounds(model)
    ocs = objective(model)
    for (i, rid) in enumerate(rxnids)
        rmets = Dict{String,Float64}()
        for (j, stoich) in zip(findnz(S[:, i])...)
            rmets[metids[j]] = stoich
        end
        modelreactions[rid] = Reaction(
            rid;
            name = reaction_name(model, rid),
            metabolites = rmets,
            lb = lbs[i],
            ub = ubs[i],
            grr = reaction_gene_association(model, rid),
            objective_coefficient = ocs[i],
            notes = reaction_notes(model, rid),
            annotations = reaction_annotations(model, rid),
            subsystem = reaction_subsystem(model, rid),
        )
    end

    return StandardModel(
        id;
        reactions = modelreactions,
        metabolites = modelmetabolites,
        genes = modelgenes,
    )
end
